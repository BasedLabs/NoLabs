import uuid
from typing import Optional, List, Generic

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import ComponentData, Job
from nolabs.domain.workflow.component import Component, ComponentTypeFactory, Parameter, TInput, TOutput
from nolabs.infrastructure.cel import cel as celery
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.red import redis_client, acquire_redlock
from nolabs.infrastructure.settings import settings
from nolabs.workflow.logic.correlation import make_correlation_id, unpack_correlation_id
from nolabs.workflow.logic.states import ControlStates, _set_state, _get_state
from nolabs.workflow.socketio_events_emitter import emit_component_jobs_event, emit_start_component_event, \
    emit_start_job_event, emit_finish_job_event, emit_finish_component_event


@celery.app.task("control._component_task", bind=True, queue=settings.workflow_queue)
async def _component_task(bind):
    correlation_id = bind.request.headers['correlation_id']
    experiment_id, component_id = unpack_correlation_id(correlation_id=correlation_id)

    emit_start_component_event(experiment_id=experiment_id, component_id=component_id)

    data: ComponentData = ComponentData.objects.with_id(component_id)

    extra = {"component_id": component_id,
             "task_id": bind.task_id,
             "correlation_id": correlation_id}

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
            correlation_id = make_correlation_id(experiment_id, component_id, job_id)

            state, _ = _get_state(id=job_id)
            if state == ControlStates.STARTED:
                continue

            redis_client.delete(correlation_id)
            task_id = uuid.UUID()
            redis_client.set(name=correlation_id, key=task_id, value=1)
            _job_task.apply_async(kwargs={
                "job_id": job_id},
                task_id=task_id,
                headers={"correlation_id": correlation_id})


async def _start_job(experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
    correlation_id = make_correlation_id(experiment_id, component_id, job_id=job_id)
    redis_client.delete(correlation_id)
    task_id = uuid.UUID()
    redis_client.set(name=correlation_id, key=task_id, value=1)
    _job_task.apply_async(kwargs={
        "job_id": job_id},
        task_id=task_id,
        headers={"correlation_id": correlation_id})


@celery.app.task("control._job_task", bind=True, queue=settings.workflow_queue)
async def _job_task(bind, job_id: uuid.UUID):
    _set_state(id=job_id, state=ControlStates.STARTED)

    correlation_id = bind.request.headers["correlation_id"]
    experiment_id, component_id = unpack_correlation_id(correlation_id=correlation_id)
    emit_start_job_event(experiment_id=experiment_id,
                         component_id=component_id,
                         job_id=job_id)
    logger.info("Job started", extra={"correlation_id": correlation_id, "job_id": job_id})
    data: ComponentData = ComponentData.objects.with_id(component_id)
    component = Component.restore(data=data)
    control_flow: ComponentFlowHandler = component.component_flow_type(
        experiment_id=experiment_id,
        component_id=component_id
    )
    await control_flow.on_job_task(job_id=job_id)


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
    async def start_component(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        correlation_id = make_correlation_id(experiment_id, component_id)

        with acquire_redlock(correlation_id):
            state, _ = _get_state(id=component_id)
            if state == ControlStates.STARTED:
                raise NoLabsException(ErrorCodes.component_running)

        await self._start_component(experiment_id=experiment_id, component_id=component_id)

    async def start_job(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
        correlation_id = make_correlation_id(experiment_id, component_id, job_id)

        with acquire_redlock(correlation_id):
            state, _ = _get_state(id=job_id)
            if state == ControlStates.STARTED:
                raise NoLabsException(ErrorCodes.job_running)

        redis_client.delete(correlation_id)
        task_id = uuid.UUID()
        redis_client.set(name=correlation_id, key=task_id, value=1)
        _job_task.apply_async(kwargs={
            "job_id": job_id},
            task_id=task_id,
            headers={"correlation_id": correlation_id})

    async def complete_job(self,
                           experiment_id: uuid.UUID,
                           component_id: uuid.UUID,
                           job_id: uuid.UUID):
        _set_state(id=job_id, state=ControlStates.SUCCESS)

        emit_finish_job_event(experiment_id=experiment_id,
                              component_id=component_id,
                              job_id=job_id)

        correlation_id = make_correlation_id(experiment_id=experiment_id,
                                             component_id=component_id,
                                             job_id=job_id)

        # acquire distributed lock at component level
        with acquire_redlock(key=correlation_id) as redlock:
            component_state, _ = _get_state(id=component_id)

            # component already succeeded, so we don't need to call completed
            if component_state == ControlStates.SUCCESS:
                return

            any_unready = False

            job_ids = Job.objects(component=component_id).only('id')

            for job_id in job_ids:
                state, _ = _get_state(id=job_id)
                if state != ControlStates.SUCCESS:
                    any_unready = True
                    break

        if not any_unready:
            await self._complete_component(experiment_id=experiment_id, component_id=component_id, job_ids=job_ids)

    async def _complete_component(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_ids: List[uuid.UUID]):
        data: ComponentData = ComponentData.objects.with_id(component_id)
        component = Component.restore(data=data)
        flow: ComponentFlowHandler = component.component_flow_type(experiment_id=experiment_id,
                                                                   component_id=component_id)
        output = await flow.on_completion(inp=component.input_value, job_ids=job_ids)
        component.output_value = output or {}
        component.dump(data=data)
        data.save()

        emit_finish_component_event(experiment_id=experiment_id, component_id=component_id)

    async def _start_component(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        correlation_id = make_correlation_id(experiment_id, component_id)

        with acquire_redlock(correlation_id):
            state, _ = _get_state(id=component_id)
            if state == ControlStates.STARTED:
                return

            _set_state(id=component_id, state=ControlStates.STARTED)

            _, matching_keys = redis_client.scan(match=f"{correlation_id}:*")

            for job_correlation_id in matching_keys:
                redis_client.delete(job_correlation_id)

            redis_client.delete(correlation_id)
            task_id = uuid.uuid4()
            redis_client.rpush(name=correlation_id, *[task_id])
            _component_task.apply_async(task_id=str(task_id),
                                        headers={'correlation_id': correlation_id})


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
        ...

    async def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
        ...

    # will be called from within task
    async def on_job_task(self, job_id: uuid.UUID):
        ...
