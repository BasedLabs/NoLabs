import asyncio
import uuid
from typing import List, Type, Any, Optional

from airflow.models import BaseOperator
from airflow.utils.context import Context
from pydantic import BaseModel

from nolabs.application.workflow import SetupOperator, OutputOperator
from nolabs.application.workflow.component import Component, JobValidationError, TOutput, TInput


class ProteinsComponentInput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinsComponentOutput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinsComponent(Component[ProteinsComponentInput, ProteinsComponentOutput]):
    @property
    def input_parameter_type(self) -> Type[TInput]:
        return ProteinsComponentInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return ProteinsComponentOutput

    @property
    def setup_operator_type(self) -> Type[BaseOperator]:
        return ProteinsSetupOperator

    @property
    def job_operator_type(self) -> Optional[Type[BaseOperator]]:
        return None

    @property
    def output_operator_type(self) -> Type[BaseOperator]:
        return ProteinsOutputOperator

    name = 'Proteins'
    description = 'Proteins datasource'


class ProteinsSetupOperator(SetupOperator):
    async def execute_async(self, context: Context) -> List[str]:
        ...


class ProteinsOutputOperator(OutputOperator):
    async def execute_async(self, context: Context) -> Any:
        component = ProteinsComponent.get(self.component_id)
        input: ProteinsComponentInput = component.input_value
        self.setup_output(ProteinsComponentOutput(proteins=input.proteins))
