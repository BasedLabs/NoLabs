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
from prefect import task, flow, State, Task, Flow
from prefect.client.schemas.objects import R, StateType, TaskRun, FlowRun
from prefect.context import get_run_context
from prefect.states import Failed

from nolabs.application.workflow.component import Component, TOutput, TInput, ComponentTypeFactory, Parameter
from nolabs.application.workflow.data import ComponentData, JobRunData
from nolabs.application.workflow.exceptions import WorkflowException, ErrorCodes


def _name_builder(name: str, id: uuid.UUID):
    return f'Name:{name},Id:{str(id)}'


def _run_name_builder(name: str, at: datetime.datetime):
    return f'{name},At:{at.isoformat()}'


class ComponentFlow(ABC, Generic[TInput, TOutput]):
    component_id: uuid.UUID
    component_name: str

    logger: logging.Logger
    '''Whether input was changed after last execution'''

    job_timeout_seconds: Optional[int] = 1
    component_timeout_seconds: Optional[int] = 10

    def __init__(self,
                 component: Component,
                 extra: Optional[Dict[str, Any]] = None):
        self.component_id = component.id
        self.component_name = component.name

        execute_flow_name = _name_builder(id=self.component_id, name=self.component_name)

        if not extra:
            self.extra = {}
        else:
            self.extra = extra

        self.execute = flow(
            name=execute_flow_name,
            flow_run_name=_run_name_builder(name=execute_flow_name, at=datetime.datetime.utcnow()),
            timeout_seconds=self.component_timeout_seconds,
            on_running=[self._on_component_running],
            on_failure=[self._on_component_completion],
            on_crashed=[self._on_component_completion],
            on_cancellation=[self._on_component_completion],
            on_completion=[self._on_component_completion]
        )(self.execute)

        execute_task_name = _name_builder(id=self.component_id, name=f'{self.component_name}-job')

        execute_task_run_name = execute_task_name + ',JobId:{job_id},At:{at}'

        self._job_task = task(
            name=execute_task_name,
            task_run_name=execute_task_run_name,
            timeout_seconds=self.job_timeout_seconds,
            on_failure=[self._on_job_completion],
            on_completion=[self._on_job_completion]
        )(self._job_task)

    async def get_jobs(self, inp: TInput) -> List[uuid.UUID]:
        return []

    async def post_execute(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
        pass

    async def _job_task(self, job_id: uuid.UUID, at: datetime.datetime) -> State[R]:
        ctx = get_run_context()
        task_run_id = ctx.task_run.id

        run = JobRunData.create(
            component_id=self.component_id,
            id=job_id,
            task_run_id=task_run_id,
            timeout=self.job_timeout_seconds,
            state=StateType.RUNNING,
            executed_at=at)
        run.save()
        return await self.job_task(job_id=job_id)

    async def job_task(self, job_id: uuid.UUID) -> State[R]:
        pass

    async def execute(self):
        data: ComponentData = ComponentData.objects.with_id(self.component_id)
        component = Component.restore(data=data)

        try:
            prev_components: List[Component] = []

            for previous_component_id in data.previous_component_ids:
                previous_component_state = ComponentData.objects.with_id(previous_component_id)
                previous_component = self._get_component(from_state=previous_component_state)

                errors = previous_component.output_errors()
                if errors:
                    return

                prev_components.append(previous_component)

            input_changed = component.set_input_from_previous(prev_components)

            if input_changed:
                input_errors = component.input_errors()

                if input_errors:
                    return

            input_value = component.input_value

            if input_changed:
                job_ids = await self.get_jobs(inp=input_value)
            else:
                job_ids = [j.id for j in JobRunData.objects(component=self.component_id).only('id')]

            job_errors = False

            JobRunData.objects(component=self.component_id).delete()

            if job_ids:
                states = self._job_task.map(job_ids, at=datetime.datetime.utcnow(), return_state=True)
                job_errors = any([s for s in states if s is Failed])

            if not job_errors:
                output = await self.post_execute(inp=input_value, job_ids=job_ids)
                component.output_value = output or {}
        except Exception as e:
            data.exception = str(e)
            raise e
        finally:
            component.dump(data=data)
            data.save()

    def _get_component(self, from_state: ComponentData) -> Component:
        component = ComponentTypeFactory.get_type(from_state.name)(
            id=from_state.id,
            input_schema=Parameter(**from_state.input_schema),
            output_schema=Parameter(**from_state.output_schema),
            input_value_dict=from_state.input_value_dict,
            output_value_dict=from_state.output_value_dict,
            previous_component_ids=from_state.previous_component_ids
        )

        return component

    def _on_component_running(self, _, flow_run: FlowRun, state: State):
        data: ComponentData = ComponentData.objects.with_id(self.component_id)
        data.flow_run_id = flow_run.id
        data.last_executed_at = datetime.datetime.utcnow()
        data.state = state.type
        data.save()

    def _on_component_completion(self, _, flow_run: FlowRun, state: State):
        data: ComponentData = ComponentData.objects(flow_run_id=flow_run.id).first()

        if state.type == StateType.CRASHED:
            pass # Can crash jobs as well

        data.state = state.type
        data.state_message = state.message
        data.save()

    def _on_job_completion(self, _, task_run: TaskRun, state: State):
        run: JobRunData = JobRunData.objects(task_run_id=task_run.id).first()
        run.state = state.type
        run.state_message = state.message
        run.save()

    async def celery_wait_async(self, async_result: AsyncResult) -> Any:
        while not async_result.ready():
            await asyncio.sleep(0.5)
        return async_result.get()
