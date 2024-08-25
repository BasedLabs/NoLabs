__all__ = [
    'ExecuteJobOperator',
    'SetupOperator',
    'OutputOperator'
]

import nest_asyncio
from mongoengine import Document, DictField, UUIDField
from pydantic import BaseModel

nest_asyncio.apply()

import uuid
from abc import abstractmethod, ABC
from typing import Any, List, Dict, Optional
from typing import Generic

import asyncio
from airflow.models import BaseOperator
from airflow.utils.context import Context
from airflow.utils.decorators import apply_defaults

from nolabs.application.workflow.component import Component, TOutput


class Communicator(Document):
    job_id: uuid.UUID = UUIDField(primary_key=True, required=True)
    input: dict = DictField()
    output: dict = DictField()


class SetupOperator(ABC, BaseOperator):
    """
    Setups component
    Setups jobs
    Emits jobs ids
    """
    component_id: uuid.UUID
    workflow_id: uuid.UUID
    _communicator_cache: List[Communicator] = []
    input_changed: bool = False
    extra: Optional[Dict[str, Any]] = None
    '''Whether input was changed after last execution'''

    @apply_defaults
    def __init__(self,
                 workflow_id: uuid.UUID,
                 component_id: uuid.UUID,
                 task_id: str,
                 extra: Optional[Dict[str, Any]] = None,
                 **kwargs):
        super().__init__(task_id=task_id, **kwargs)

        self.component_id = component_id
        self.workflow_id = workflow_id
        self.extra = extra

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

    def post_execute(self, context: Any, result: Any = None):
        for c in self._communicator_cache:
            c.save()

    def get_component(self) -> Component:
        return Component.get(self.component_id)

    def set_job_input(self, job_id: uuid.UUID, data: BaseModel):
        doc = Communicator.objects.with_id(job_id)

        value = data.dict()

        if not doc:
            doc = Communicator(job_id=job_id, value=value)
        doc.value = value
        self._communicator_cache.append(doc)


class ExecuteJobOperator(ABC, BaseOperator):
    component_id: uuid.UUID
    job_id: uuid.UUID
    extra: Optional[Dict[str, Any]] = None,
    _communicator_cache: Optional[Communicator]

    @apply_defaults
    def __init__(self,
                 workflow_id: uuid.UUID,
                 job_id: str,
                 component_id: uuid.UUID,
                 task_id: str,
                 extra: Optional[Dict[str, Any]] = None,
                 **kwargs):
        super().__init__(task_id=task_id, **kwargs)

        self.component_id = component_id
        self.workflow_id = workflow_id
        self.job_id = uuid.UUID(job_id)

        self._refresh_communicator()
        self.extra = extra

    @abstractmethod
    async def execute_async(self, context: Context) -> Any:
        ...

    def execute(self, context: Context) -> Any:
        loop = asyncio.get_event_loop()
        return loop.run_until_complete(self.execute_async(context))

    def post_execute(self, context: Context, result: Any = None):
        if self._communicator_cache and (self._communicator_cache.input or self._communicator_cache.output):
            self._communicator_cache.save()

    def _refresh_communicator(self):
        if not self._communicator_cache:
            self._communicator_cache = Communicator(job_id=self.job_id)
        else:
            self._communicator_cache = Communicator.objects.with_id(self.job_id)

    def get_component(self) -> Component:
        return Component.get(self.component_id)

    def get_input(self) -> Optional[Dict[str, Any]]:
        return self._communicator_cache.input

    def set_output(self, value: BaseModel):
        self._communicator_cache.output = value.dict()


class OutputOperator(ABC, BaseOperator, Generic[TOutput]):
    """
    Jobs post-processing
    """
    component_id: uuid.UUID
    workflow_id: uuid.UUID
    setup_output_called: bool = False
    extra: Optional[Dict[str, Any]] = None

    @apply_defaults
    def __init__(self, workflow_id: uuid.UUID, component_id: uuid.UUID, task_id: str,
                 extra: Optional[Dict[str, Any]] = None,
                 **kwargs):
        super().__init__(task_id=task_id, **kwargs)

        self.component_id = component_id
        self.workflow_id = workflow_id
        self.extra = extra

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

    def get_component(self) -> Component:
        return Component.get(self.component_id)

    def get_job_output(self, job_id: uuid.UUID) -> Optional[Dict[str, Any]]:
        doc: Communicator = Communicator.objects.with_id(job_id)

        if not doc or not doc.output:
            return None

        return doc.output
