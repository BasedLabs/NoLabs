import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.application.workflow.component import Component, JobValidationError


class ProteinsComponentInput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinsComponentOutput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinsComponent(Component[ProteinsComponentInput, ProteinsComponentOutput]):
    name = 'Proteins'
    description = 'Proteins datasource'

    async def execute(self):
        self.output = ProteinsComponentOutput(
            proteins=self.input.proteins
        )

    async def setup_jobs(self):
        pass

    async def jobs_setup_errors(self) -> List[JobValidationError]:
        return []

    @property
    def _input_parameter_type(self) -> Type[ProteinsComponentInput]:
        return ProteinsComponentInput

    @property
    def _output_parameter_type(self) -> Type[ProteinsComponentOutput]:
        return ProteinsComponentOutput
