__all__ = [
    'ExecuteJobTask',
    'SetupTask',
    'OutputTask'
]

import asyncio
import datetime
import logging

import nest_asyncio
from celery.result import AsyncResult
from prefect import task, get_run_logger
from prefect.context import get_run_context
from pydantic import BaseModel

import uuid
from abc import abstractmethod, ABC
from typing import Any, List, Dict, Optional
from typing import Generic

from nolabs.application.workflow.component import Component, TOutput, TInput


def _name_builder(name: str, id: uuid.UUID, t: str, method: str):
    return f'Name:{name},Id:{str(id)},Task:{t},Method:{method}'


def _run_name_builder(name: str, at: datetime.datetime):
    return f'{name},At:{at.isoformat()}'


class SetupTask(ABC, Generic[TInput]):
    """
    Setups component
    Setups jobs
    Emits jobs ids
    """
    component_id: uuid.UUID
    workflow_id: uuid.UUID
    input_changed: bool = False
    extra: Optional[Dict[str, Any]] = None

    logger: logging.Logger
    '''Whether input was changed after last execution'''

    timeout_seconds: Optional[float] = 1.0

    def __init__(self,
                 workflow_id: uuid.UUID,
                 component_id: uuid.UUID,
                 component_name: str,
                 extra: Optional[Dict[str, Any]] = None):
        self.component_id = component_id
        self.component_name = component_name
        self.workflow_id = workflow_id
        self.extra = extra
        self.logger = get_run_logger()

        setup_task_name = _name_builder(name=component_name, id=component_id, t='Setup',
                                        method=self.execute.__name__)

        self.execute_task = task(
            name=setup_task_name,
            task_run_name=_run_name_builder(name=setup_task_name, at=datetime.datetime.utcnow()),
            timeout_seconds=self.timeout_seconds
        )(self._execute)

    @abstractmethod
    async def execute(self) -> List[uuid.UUID]:
        ...

    async def _execute(self) -> List[uuid.UUID]:
        ctx = get_run_context()
        task_run_id = ctx.task_run.id

        component = self.get_component()

        prev_components: List[Component] = []

        for previous_component_id in component.previous_component_ids:
            previous_component = Component.get(previous_component_id)

            errors = previous_component.output_errors()
            if errors:
                raise ValueError(errors[0].msg)

            prev_components.append(previous_component)

        self.input_changed = component.set_input_from_previous(prev_components)

        errors = component.input_errors()
        if errors:
            raise ValueError(errors[0].msg)

        component.save()

        return await self.execute()

    def get_component(self) -> Component:
        return Component.get(self.component_id)


class ExecuteJobTask(ABC):
    component_id: uuid.UUID
    extra: Optional[Dict[str, Any]] = None,
    logger: logging.Logger
    timeout_seconds: Optional[float] = 10.0

    def __init__(self,
                 workflow_id: uuid.UUID,
                 component_name: str,
                 component_id: uuid.UUID,
                 extra: Optional[Dict[str, Any]] = None):
        self.component_id = component_id
        self.workflow_id = workflow_id
        self.logger = get_run_logger()

        execute_task_name = _name_builder(name=component_name, id=component_id, t='ExecuteJob',
                                          method=self.execute.__name__)

        self.execute_task = task(
            name=execute_task_name,
            task_run_name=_run_name_builder(name=execute_task_name, at=datetime.datetime.utcnow()),
            timeout_seconds=self.timeout_seconds
        )(self._execute)

        self.extra = extra

    @abstractmethod
    async def execute(self, job_id: uuid.UUID) -> Optional[BaseModel]:
        ...

    async def _execute(self, job_id: uuid.UUID):
        await self.execute(job_id=job_id)

    def get_component(self) -> Component:
        return Component.get(self.component_id)

    async def celery_wait_async(self, async_result: AsyncResult) -> Any:
        while not async_result.ready():
            await asyncio.sleep(0.5)
        return async_result.get()


class OutputTask(ABC, Generic[TOutput]):
    """
    Jobs post-processing
    """
    workflow_id: uuid.UUID
    component_id: uuid.UUID
    setup_output_called = False
    extra: Optional[Dict[str, Any]] = None
    timeout_seconds: Optional[float] = 1.0

    def __init__(self, workflow_id: uuid.UUID,
                 component_name: str,
                 component_id: uuid.UUID,
                 extra: Optional[Dict[str, Any]] = None):
        self.workflow_id = workflow_id
        self.component_id = component_id

        execute_task_name = _name_builder(name=component_name, id=component_id, t='Output',
                                          method=self.execute.__name__)

        self.execute_task = task(
            name=execute_task_name,
            task_run_name=_run_name_builder(name=execute_task_name, at=datetime.datetime.utcnow()),
            timeout_seconds=self.timeout_seconds
        )(self.execute)

        post_execute_name = _name_builder(name=component_name, id=component_id, t='Output',
                                          method=self.post_execute.__name__)

        self.post_execute_task = task(
            name=post_execute_name,
            task_run_name=_run_name_builder(name=post_execute_name, at=datetime.datetime.utcnow()),
            timeout_seconds=self.timeout_seconds
        )(self.post_execute)

        self.extra = extra

    @abstractmethod
    async def execute(self) -> Optional[BaseModel]:
        ...

    async def post_execute(self):
        if not self.setup_output_called:
            raise ValueError('Setup output not called')

    def setup_output(self, output: TOutput):
        component = Component.get(self.component_id)
        component.output_value = output
        component.save()
        self.setup_output_called = True

    def get_component(self) -> Component:
        return Component.get(self.component_id)
