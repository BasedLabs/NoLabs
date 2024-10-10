import asyncio
import uuid
from typing import Optional, List, Generic, Dict, Any

from nolabs.domain.models.common import ComponentData
from nolabs.infrastructure.cel import cel as celery
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.settings import settings
from nolabs.workflow.component import Component, ComponentTypeFactory, Parameter, TInput, TOutput


@celery.app.task(name="control._component_main_task", bind=True, queue=settings.workflow_queue)
def _component_main_task(bind, experiment_id: uuid.UUID, component_id: uuid.UUID) -> List[uuid.UUID]:
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
        job_ids = await control_flow.on_component_task(inp=input_value)

        logger.info("Retrieved job ids", extra={**extra, **{"job_ids": job_ids}})

        return job_ids

    asyncio.run(_())


@celery.app.task(name="control._job_main_task", bind=True, queue=settings.workflow_queue)
def _job_main_task(bind, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
    async def _():
        data: ComponentData = ComponentData.objects.with_id(component_id)
        component = Component.restore(data=data)
        control_flow: ComponentFlowHandler = component.component_flow_type(
            experiment_id=experiment_id,
            component_id=component_id
        )
        await control_flow.on_job_task(job_id=job_id)

    asyncio.run(_())


@celery.app.task(name="control._complete_job_task", bind=True, queue=settings.workflow_queue)
def _complete_job_task(bind,
                       experiment_id: uuid.UUID,
                       component_id: uuid.UUID,
                       job_id: uuid.UUID,
                       long_running_output: Optional[Dict[str, Any]]):
    async def _():
        data: ComponentData = ComponentData.objects.with_id(component_id)
        component = Component.restore(data=data)
        control_flow: ComponentFlowHandler = component.component_flow_type(
            experiment_id=experiment_id,
            component_id=component_id
        )
        await control_flow.on_job_task(job_id=job_id, long_running_output=long_running_output)

    asyncio.run(_())


@celery.app.task(name="control._complete_component_task", bind=True, queue=settings.workflow_queue)
async def _complete_component_task(experiment_id: uuid.UUID, component_id: uuid.UUID, job_ids: List[uuid.UUID]):
    data: ComponentData = ComponentData.objects.with_id(component_id)
    component = Component.restore(data=data)
    flow: ComponentFlowHandler = component.component_flow_type(experiment_id=experiment_id,
                                                               component_id=component_id)
    output = await flow.on_completion(inp=component.input_value, job_ids=job_ids)
    component.output_value = output or {}
    component.dump(data=data)
    data.save()

    logger.info("Component finished", extra={"component_id": component_id})


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


class ComponentFlowHandler(Generic[TInput, TOutput]):
    """
    Represents component flow and available client handlers
    """

    def __init__(self,
                 component_id: uuid.UUID,
                 experiment_id: uuid.UUID):
        self.component_id = component_id
        self.experiment_id = experiment_id

    async def on_component_task(self, inp: TInput) -> List[uuid.UUID]:
        return []

    async def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
        return None

    async def on_job_task(self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]):
        pass

    async def on_job_completion(self, job_id: uuid.UUID):
        pass
