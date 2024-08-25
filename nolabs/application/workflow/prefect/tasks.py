__all__ = [
    'ExecuteJobTask',
    'SetupTask',
    'OutputTask'
]

import logging

import nest_asyncio
import pydantic
from prefect import task, get_run_logger
from pydantic import BaseModel

nest_asyncio.apply()

import uuid
from abc import abstractmethod, ABC
from typing import Any, List, Dict, Optional, TYPE_CHECKING
from typing import Generic

from nolabs.application.workflow.component import Component, TOutput, TInput


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

    def __init__(self,
                 workflow_id: uuid.UUID,
                 component_id: uuid.UUID,
                 extra: Optional[Dict[str, Any]] = None):
        self.component_id = component_id
        self.workflow_id = workflow_id
        self.extra = extra
        self.logger = get_run_logger()

    @task(description='Setups flow')
    def pre_execute(self):
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

    @abstractmethod
    async def _execute(self) -> List[uuid.UUID]:
        """
                Setups jobs
                Returns list of job ids
        """
        ...

    @task(description='Main setup task execution', name='{component_name}', task_run_name='{component_id}')
    async def execute(self, component_name: str, component_id: uuid.UUID) -> List[uuid.UUID]:
        self.pre_execute.submit()
        return await self._execute()

    def get_component(self) -> Component:
        return Component.get(self.component_id)


class ExecuteJobTask(ABC):
    component_id: uuid.UUID
    extra: Optional[Dict[str, Any]] = None,
    logger: logging.Logger

    def __init__(self,
                 workflow_id: uuid.UUID,
                 component_id: uuid.UUID,
                 extra: Optional[Dict[str, Any]] = None):
        self.component_id = component_id
        self.workflow_id = workflow_id
        self.logger = get_run_logger()
        self.extra = extra

    @abstractmethod
    async def _execute(self, job_id: uuid.UUID) -> Optional[BaseModel]:
        ...

    @task(description='Executes a job', name='{component_name}', task_run_name='{component_id}')
    async def execute(self, job_id: uuid.UUID, component_name: str, component_id: uuid.UUID) -> Optional[BaseModel]:
        return await self._execute(job_id=job_id)

    def get_component(self) -> Component:
        return Component.get(self.component_id)


class OutputTask(ABC, Generic[TOutput]):
    """
    Jobs post-processing
    """
    workflow_id: uuid.UUID
    component_id: uuid.UUID
    setup_output_called = False
    extra: Optional[Dict[str, Any]] = None

    def __init__(self, workflow_id: uuid.UUID, component_id: uuid.UUID, extra: Optional[Dict[str, Any]] = None):
        self.workflow_id = workflow_id
        self.component_id = component_id
        self.extra = extra

    @abstractmethod
    async def _execute(self) -> Optional[BaseModel]:
        """
            Post jobs processing
            Setup component output data
        """
        ...

    @task(description='Setups output', name='{component_name}', task_run_name='{component_id}')
    async def execute(self, component_name: str, component_id: uuid.UUID) -> Optional[BaseModel]:
        result = await self._execute()

        if not self.setup_output_called:
            raise ValueError('Setup output not called')

        return result

    def setup_output(self, output: TOutput):
        component = Component.get(self.component_id)
        component.output_value = output
        component.save()
        self.setup_output_called = True

    def get_component(self) -> Component:
        return Component.get(self.component_id)
