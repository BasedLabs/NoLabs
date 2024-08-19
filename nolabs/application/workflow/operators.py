__all__ = [
    'JobOperator',
    'SetupOperator',
    'OutputOperator'
]

import nest_asyncio
nest_asyncio.apply()

import uuid
from abc import abstractmethod, ABC
from typing import Any, List
from typing import Generic

import asyncio
from airflow.models import BaseOperator
from airflow.utils.context import Context
from airflow.utils.decorators import apply_defaults

from nolabs.application.workflow.component import Component, TOutput
from nolabs.domain.models.common import JobId



class SetupOperator(ABC, BaseOperator):
    """
    Setups component
    Setups jobs
    Emits jobs ids
    """
    component_id: uuid.UUID
    workflow_id: uuid.UUID
    experiment_id: uuid.UUID
    input_changed: bool = False
    '''Whether input was changed after last execution'''

    @apply_defaults
    def __init__(self, workflow_id: uuid.UUID, experiment_id: uuid.UUID, component_id: uuid.UUID, task_id: str,
                 **kwargs):
        super().__init__(task_id=task_id, **kwargs)

        self.component_id = component_id
        self.workflow_id = workflow_id
        self.experiment_id = experiment_id

    def pre_execute(self, context: Any):
        component = Component.get(self.component_id)

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
    async def execute_async(self, context: Context) -> List[str]:
        """
                Setups jobs
                Returns list of job ids
        """
        ...

    def execute(self, context: Context) -> List[str]:
        loop = asyncio.get_event_loop()
        return loop.run_until_complete(self.execute_async(context))

    def serialize_job_ids(self, job_ids: List[JobId]) -> List[str]:
        return [str(i) for i in job_ids]


class JobOperator(ABC, BaseOperator):
    """
    Executes job
    """
    component_id: uuid.UUID
    job_id: JobId

    @apply_defaults
    def __init__(self, workflow_id: uuid.UUID, experiment_id: uuid.UUID, job_id: str, component_id: uuid.UUID,
                 task_id: str, **kwargs):
        super().__init__(task_id=task_id, **kwargs)

        self.component_id = component_id
        self.workflow_id = workflow_id
        self.experiment_id = experiment_id
        self.job_id = JobId(uuid.UUID(job_id))

    @abstractmethod
    async def execute_async(self, context: Context) -> Any:
        ...

    def execute(self, context: Context) -> Any:
        loop = asyncio.get_event_loop()
        return loop.run_until_complete(self.execute_async(context))


class OutputOperator(ABC, BaseOperator, Generic[TOutput]):
    """
    Jobs post-processing
    """
    component_id: uuid.UUID
    workflow_id: uuid.UUID
    experiment_id: uuid.UUID
    setup_output_called: bool = False

    @apply_defaults
    def __init__(self, workflow_id: uuid.UUID, experiment_id: uuid.UUID, component_id: uuid.UUID, task_id: str,
                 **kwargs):
        super().__init__(task_id=task_id, **kwargs)

        self.component_id = component_id
        self.workflow_id = workflow_id
        self.experiment_id = experiment_id

    @abstractmethod
    async def execute_async(self, context: Context) -> Any:
        """
            Post jobs processing
            Setup component output data
        """
        ...

    def execute(self, context: Context) -> Any:
        loop = asyncio.get_event_loop()
        result = loop.run_until_complete(self.execute_async(context))

        if not self.setup_output_called:
            raise ValueError('You must setup output')

        return result

    def setup_output(self, output: TOutput):
        component = Component.get(self.component_id)
        component.output_value = output
        component.save()
        self.setup_output_called = True
