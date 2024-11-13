import asyncio
import uuid
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from asgiref.sync import async_to_sync
from celery import Celery
from celery.result import allow_join_result

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import ComponentData, Experiment, Job
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.redis_client_factory import redlock
from nolabs.infrastructure.settings import settings
from nolabs.workflow.core import Tasks
from nolabs.workflow.core.component import Component, ComponentTypeFactory, Parameter
from nolabs.workflow.core.graph import Graph

if TYPE_CHECKING:
    from nolabs.workflow.core.flow import ComponentFlowHandler


def register_workflow_celery_tasks(celery: Celery):
    @celery.task(
        name=Tasks.component_main_task,
        bind=True,
        queue=settings.workflow_queue,
        max_retries=0,
    )
    def component_main_task(bind, experiment_id: uuid.UUID, component_id: uuid.UUID):
        async def _():
            data: ComponentData = ComponentData.objects.with_id(component_id)

            extra = {"component_id": component_id, "celery_task_id": bind.request.id}

            component = Component.restore(data=data)
            prev_components: List[Component] = []

            for previous_component_id in data.previous_component_ids:
                previous_component_state = ComponentData.objects.with_id(
                    previous_component_id
                )
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

            component.set_input_from_previous(prev_components)
            input_errors = component.input_errors()

            if input_errors:
                logger.info(
                    "Input errors",
                    extra={
                        **extra,
                        **{"input_errors": [(e.msg, e.loc) for e in input_errors]},
                    },
                )
                raise NoLabsException(ErrorCodes.component_input_invalid, 'Invalid component input')

            component.dump(data=data)
            data.save()

            control_flow: ComponentFlowHandler = component.component_flow_type(
                experiment_id=experiment_id, component_id=component_id
            )

            input_value = component.input_value
            await control_flow.on_component_task(inp=input_value)

        async_to_sync(_)()

    @celery.task(
        name="workflow._job_main_task",
        bind=True,
        queue=settings.workflow_queue,
        max_retries=0,
    )
    def job_main_task(
        bind, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID
    ):
        async def _():
            data: ComponentData = ComponentData.objects.with_id(component_id)
            component = Component.restore(data=data)
            control_flow: ComponentFlowHandler = component.component_flow_type(
                experiment_id=experiment_id, component_id=component_id
            )
            with allow_join_result():
                await control_flow.on_job_task(job_id=job_id)

        async_to_sync(_)()

    @celery.task(
        name=Tasks.complete_job_task,
        bind=True,
        queue=settings.workflow_queue,
        max_retries=0,
    )
    def complete_job_task(
        bind,
        experiment_id: uuid.UUID,
        component_id: uuid.UUID,
        job_id: uuid.UUID,
        long_running_output: Optional[Dict[str, Any]],
    ):
        async def _():
            data: ComponentData = ComponentData.objects.with_id(component_id)
            component = Component.restore(data=data)
            control_flow: ComponentFlowHandler = component.component_flow_type(
                experiment_id=experiment_id, component_id=component_id
            )
            await control_flow.on_job_completion(
                job_id=job_id, long_running_output=long_running_output
            )

        asyncio.run(_())

    @celery.task(
        name=Tasks.complete_component_task,
        bind=True,
        queue=settings.workflow_queue,
        max_retries=0,
    )
    def complete_component_task(
        bind, experiment_id: uuid.UUID, component_id: uuid.UUID
    ):
        async def _():
            data: ComponentData = ComponentData.objects.with_id(component_id)
            component = Component.restore(data=data)
            flow: ComponentFlowHandler = component.component_flow_type(
                experiment_id=experiment_id, component_id=component_id
            )
            jobs = Job.objects(component=component_id).only("id")
            output = await flow.on_completion(
                inp=component.input_value, job_ids=[j.id for j in jobs]
            )
            component.output_value = output or {}
            component.dump(data=data)
            data.save()

            logger.info("Component finished", extra={"component_id": component_id})

        async_to_sync(_)()

    @celery.task(name=Tasks.sync_graphs_task, bind=True, queue=settings.workflow_queue)
    def sync_graphs(bind):
        async def _():
            lock = redlock(key=Tasks.sync_graphs_task, auto_release_time=10.0)

            if not lock.acquire():
                return

            try:
                experiments = Experiment.objects.only("id")

                for experiment in experiments:
                    graph = Graph(experiment_id=experiment.id)
                    await graph.sync()
            finally:
                if lock.locked():
                    lock.release()

        async_to_sync(_)()


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
