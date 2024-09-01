__all__ = [
    'ComponentFlow'
]

import asyncio
import datetime
import logging
import uuid
from abc import ABC
from typing import Any, List, Dict, Optional
from typing import Generic

from celery.result import AsyncResult
from prefect import task, flow
from prefect.context import get_run_context

from nolabs.application.workflow.component import Component, TOutput, TInput, ComponentTypeFactory, Parameter
from nolabs.application.workflow.data import ComponentState, ComponentRunModel, JobRunModel, InputPropertyErrorModel


def _name_builder(name: str, id: uuid.UUID):
    return f'Name:{name},Id:{str(id)}'


def _run_name_builder(name: str, at: datetime.datetime):
    return f'{name},At:{at.isoformat()}'


class ComponentFlow(ABC, Generic[TInput, TOutput]):
    component_id: uuid.UUID
    workflow_id: uuid.UUID
    input_changed: bool = False
    force_run: bool = False
    extra: Optional[Dict[str, Any]] = None

    logger: logging.Logger
    '''Whether input was changed after last execution'''

    job_timeout_seconds: Optional[float] = 1.0
    component_timeout_seconds: Optional[float] = 10.0

    def __init__(self,
                 component_id: uuid.UUID,
                 component_name: str,
                 workflow_id: uuid.UUID,
                 force_run: bool = False,
                 extra: Optional[Dict[str, Any]] = None):
        self.component_id = component_id
        self.workflow_id = workflow_id
        self.force_run = force_run

        execute_flow_name = _name_builder(id=component_id, name=component_name)

        if not extra:
            self.extra = {}
        else:
            self.extra = extra

        self.execute = flow(
            name=execute_flow_name,
            flow_run_name=_run_name_builder(name=execute_flow_name, at=datetime.datetime.utcnow()),
            timeout_seconds=self.component_timeout_seconds
        )(self.execute)

        execute_task_name = _name_builder(id=component_id, name=f'{component_name}-job')

        execute_task_run_name = execute_task_name + ',JobId:{job_id},At:{at}'

        self._job_task = task(
            name=execute_task_name,
            task_run_name=execute_task_run_name,
            timeout_seconds=self.job_timeout_seconds
        )(self._job_task)

    async def pre_execute(self, inp: TInput):
        pass

    async def get_jobs(self, inp: TInput) -> List[uuid.UUID]:
        return []

    async def post_execute(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
        pass

    async def _job_task(self, job_id: uuid.UUID, at: datetime.datetime):
        state: ComponentState = ComponentState.objects.with_id(self.component_id)
        await self._add_job_task_run(job_id=job_id, state=state)
        state.save()
        await self.job_task(job_id=job_id)

    async def job_task(self, job_id: uuid.UUID):
        pass

    async def execute(self):
        state: ComponentState = ComponentState.objects.with_id(self.component_id)

        self._add_task_run(state)

        state.save()

        component = self._get_component(from_state=state)

        prev_components: List[Component] = []

        for previous_component_id in state.previous_component_ids:
            previous_component_state = ComponentState.objects.with_id(previous_component_id)
            previous_component = self._get_component(from_state=previous_component_state)

            errors = previous_component.output_errors()
            if errors:
                raise ValueError('There were errors in previous component outputs')

            prev_components.append(previous_component)

        input_changed = component.set_input_from_previous(prev_components)

        if not input_changed and not self.force_run:
            return

        errors = component.input_errors()
        if errors:
            state.input_property_errors = [InputPropertyErrorModel.create(loc=err.loc, msg=err.msg) for err in errors]
            state.save()
            raise ValueError('There were errors in this component inputs')

        state.set_component(component=component)
        state.save()

        input_value = component.input_value

        await self.pre_execute(input_value)

        job_ids = await self.get_jobs(inp=input_value)

        if job_ids:
            self._job_task.map(job_ids, at=datetime.datetime.utcnow()).wait()

        output = await self.post_execute(inp=input_value, job_ids=job_ids)
        if output:
            component.output_value = output
            state.set_component(component=component)
            state.save()

    def _get_component(self, from_state: ComponentState) -> Component:
        component = ComponentTypeFactory.get_type(from_state.name)(
            id=from_state.id,
            input_schema=Parameter(**from_state.input_schema),
            output_schema=Parameter(**from_state.output_schema),
            input_value_dict=from_state.input_value_dict,
            output_value_dict=from_state.output_value_dict,
            previous_component_ids=from_state.previous_component_ids
        )

        return component

    @staticmethod
    async def celery_wait_async(async_result: AsyncResult) -> Any:
        while not async_result.ready():
            await asyncio.sleep(0.5)
        return async_result.get()

    def _add_task_run(self, state: ComponentState):
        ctx = get_run_context()
        task_run_id = ctx.flow_run.id

        previous_run = state.runs[-2] if len(state.runs) > 1 else None

        state.runs.append(
            ComponentRunModel.create(task_run_id=task_run_id,
                                     jobs=previous_run.jobs if previous_run else [],
                                     created_at=datetime.datetime.utcnow()))

    def _add_job_task_run(self, job_id: uuid.UUID, state: ComponentState):
        ctx = get_run_context()
        task_run_id = ctx.task_run.id

        if not state.runs:
            raise ValueError('Runs were not found in component')

        run = state.runs[-1]

        run.jobs.append(JobRunModel.create(id=job_id, task_run_id=task_run_id, created_at=datetime.datetime.utcnow()))
