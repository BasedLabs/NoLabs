import uuid
from typing import List, Type, Any, Optional

from prefect import task
from pydantic import BaseModel

from nolabs.application.workflow import Component, SetupTask, OutputTask, ExecuteJobTask
from nolabs.application.workflow.component import TInput, TOutput


class ProteinsComponentInput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinsComponentOutput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinsComponent(Component[ProteinsComponentInput, ProteinsComponentOutput]):
    name = 'Proteins'
    description = 'Proteins datasource'

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return ProteinsComponentInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return ProteinsComponentOutput

    @property
    def setup_task_type(self) -> Type[SetupTask]:
        return ProteinsSetupTask

    @property
    def job_task_type(self) -> Optional[Type[ExecuteJobTask]]:
        return None

    @property
    def output_task_type(self) -> Type[OutputTask]:
        return ProteinsOutputTask


class ProteinsSetupTask(SetupTask):

    async def execute(self) -> List[uuid.UUID]:
        return []


class ProteinsOutputTask(OutputTask):
    timeout_seconds = 1.0

    async def execute(self) -> Optional[BaseModel]:
        component = ProteinsComponent.get(self.component_id)
        input: ProteinsComponentInput = component.input_value
        self.setup_output(ProteinsComponentOutput(proteins=input.proteins))
        return
