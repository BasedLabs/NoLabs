import asyncio
import uuid
from typing import Optional, List, Generic

from celery.result import AsyncResult

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import ComponentData
from nolabs.infrastructure.cel import cel as celery
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.red import redis_client, redlock
from nolabs.infrastructure.settings import settings
from nolabs.workflow.component import Component, ComponentTypeFactory, Parameter, TInput, TOutput
from nolabs.workflow.logic.correlation import _assign_correlation_id, _make_component_cid, _make_job_cid
from nolabs.workflow.logic.states import _ready, ControlStates, _set_internal_state
from nolabs.workflow.socketio_events_emitter import emit_component_jobs_event, emit_start_component_event, \
    emit_finish_job_event


@celery.app.task(name="control._component_task", bind=True, queue=settings.workflow_queue)
def _component_task(bind, experiment_id: uuid.UUID, component_id: uuid.UUID):
    async def _():
        data: ComponentData = ComponentData.objects.with_id(component_id)

        extra = {"component_id": component_id,
                 "celery_task_id": bind.request.id}

        component = Component.restore(data=data)
        prev_components: List[Component] = []

        for previous_component_id in data.previous_component_ids:
            previous_component_state = ComponentData.objects.with_id(previous_component_id)
            previous_component = _get_component(from_state=previous_component_state)

            errors = previous_component.output_errors()
            if errors:
                logger.info(
                    "Previous component output errors",
                    extra={
                        **extra,
                        **{
                            "previous_component_id": previous_component_id,
                            "errors": [(e.msg, e.loc) for e in errors],
                        },
                    },
                )
                return

            prev_components.append(previous_component)

        input_changed = component.set_input_from_previous(prev_components)

        if input_changed:
            logger.info("Input changed, checking input", extra=extra)

            input_errors = component.input_errors()

            if input_errors:
                logger.info(
                    "Input errors",
                    extra={
                        **extra,
                        **{"input_errors": [(e.msg, e.loc) for e in input_errors]},
                    },
                )
                return

        control_flow: ComponentFlowHandler = component.component_flow_type(
            experiment_id=experiment_id,
            component_id=component_id
        )

        input_value = component.input_value
        job_ids = await control_flow.on_started(inp=input_value)

        logger.info("Retrieved job ids", extra={**extra, **{"job_ids": job_ids}})

        if job_ids:
            emit_component_jobs_event(
                experiment_id=experiment_id,
                component_id=component_id,
                job_ids=job_ids
            )

            for job_id in job_ids:
                await _start_job(job_id=job_id)

    asyncio.run(_())


async def _start_job(job_id: uuid.UUID):
    lock = redlock(key=str(job_id))
    try:
        # Locking with blocking=False to prevent simultaneous running
        # Skip execution if someone is executing it right now
        if lock.acquire(blocking=False):
            if not await _ready(id=job_id):
                logger.info("Cannot start job since it is already executing",
                            extra={"job_id": job_id})
                return

            correlation_id = _make_job_cid(job_id)
            celery_task_id = await _assign_correlation_id(correlation_id)
            _job_task.apply_async(kwargs={
                "job_id": job_id},
                task_id=celery_task_id,
                retry=False)
    finally:
        lock.release()


@celery.app.task(name="control._job_task", bind=True, queue=settings.workflow_queue)
def _job_task(bind, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
    async def _():
        data: ComponentData = ComponentData.objects.with_id(component_id)
        component = Component.restore(data=data)
        control_flow: ComponentFlowHandler = component.component_flow_type(
            experiment_id=experiment_id,
            component_id=component_id
        )
        await control_flow.on_job_task(job_id=job_id)

    asyncio.run(_())


async def _complete_component(experiment_id: uuid.UUID, component_id: uuid.UUID, job_ids: List[uuid.UUID]):
    data: ComponentData = ComponentData.objects.with_id(component_id)
    component = Component.restore(data=data)
    flow: ComponentFlowHandler = component.component_flow_type(experiment_id=experiment_id,
                                                               component_id=component_id)
    output = await flow.on_completion(inp=component.input_value, job_ids=job_ids)
    component.output_value = output or {}
    component.dump(data=data)
    data.save()

    logger.info("Component finished", extra={"component_id": component_id})

    #next_components = ComponentData.objects(previous_component_ids=component_id).only('id')
#
    #for next_component in next_components:
    #    await _start_component(experiment_id=experiment_id, component_id=next_component.id)


async def _start_component(experiment_id: uuid.UUID, component_id: uuid.UUID):
    correlation_id = _make_component_cid(component_id)

    lock = redlock(key=correlation_id)
    try:
        # Locking with blocking=False to prevent simultaneous running
        # Skip execution if someone is executing it right now
        if lock.acquire(blocking=False):
            if not await _ready(id=component_id):
                logger.info("Cannot start component since it is already executing",
                            extra={"correlation_id": correlation_id, "component_id": component_id})
                return

            _, matching_keys = await redis_client.scan(match=f"{correlation_id}*")

            for job_correlation_id in matching_keys:
                await redis_client.delete(job_correlation_id)

            celery_task_id = await _assign_correlation_id(correlation_id=correlation_id)
            _component_task.apply_async(task_id=celery_task_id, kwargs={
                "experiment_id": experiment_id,
                "component_id": component_id
            }, retry=False)
            emit_start_component_event(experiment_id=experiment_id,
                                       component_id=component_id)
            logger.info("Component started", extra={"component_id": component_id})
    finally:
        lock.release()


def _get_component(from_state: ComponentData) -> Component:
    component = ComponentTypeFactory.get_type(from_state.name)(
        id=from_state.id,
        input_schema=Parameter(**from_state.input_schema),
        output_schema=Parameter(**from_state.output_schema),
        input_value_dict=from_state.input_value_dict,
        output_value_dict=from_state.output_value_dict,
        previous_component_ids=from_state.previous_component_ids,
    )

    return component


class Flow:
    @classmethod
    async def start_component(cls, experiment_id: uuid.UUID, component_id: uuid.UUID):

        lock = redlock(key=str(component_id))
        try:
            if lock.acquire(blocking=False):
                if not await _ready(id=component_id):
                    raise NoLabsException(ErrorCodes.component_running)

                await _start_component(experiment_id=experiment_id, component_id=component_id)
        finally:
            lock.release()

    @classmethod
    async def start_job(cls, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
        lock = redlock(key=str(job_id))
        try:
            if lock.acquire(blocking=False):
                if not await _ready(id=component_id):
                    raise NoLabsException(ErrorCodes.job_running)

                await _start_job(job_id=job_id)
        finally:
            lock.release()


class ComponentFlowHandler(Generic[TInput, TOutput]):
    """
    Represents component flow and available client handlers
    """

    def __init__(self,
                 component_id: uuid.UUID,
                 experiment_id: uuid.UUID):
        self.component_id = component_id
        self.experiment_id = experiment_id

    # will be called from within task
    async def on_started(self, inp: TInput) -> List[uuid.UUID]:
        return []

    async def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
        return None

    # will be called from within task
    async def on_job_task(self, job_id: uuid.UUID):
        pass

    async def complete_job(self, job_id: uuid.UUID):
        return await self._complete_job(job_id=job_id, state=ControlStates.SUCCESS)

    async def cancel_job(self, job_id: uuid.UUID):
        return await self._complete_job(job_id=job_id, state=ControlStates.CANCELLED)

    async def _complete_job(self, job_id: uuid.UUID, state: ControlStates, **kwargs):
        await _set_internal_state(id=job_id, state=state, state_message=kwargs.get("state_message"))

        emit_finish_job_event(experiment_id=self.experiment_id,
                              component_id=self.component_id,
                              job_id=job_id)
