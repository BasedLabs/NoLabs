import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.workflow.component import PythonComponent, JobValidationError


class ProteinsComponentInput(BaseModel):
    protein_ids: List[uuid.UUID]


class ProteinsComponentOutput(BaseModel):
    protein_ids: List[uuid.UUID]


class ProteinsComponent(PythonComponent[ProteinsComponentInput, ProteinsComponentOutput]):
    name = 'Proteins'

    async def execute(self):
        self.output = ProteinsComponentOutput(
            protein_ids=self.input.protein_ids
        )

    async def setup_jobs(self):
        pass

    async def prevalidate_jobs(self) -> List[JobValidationError]:
        return []

    @property
    def _input_parameter_type(self) -> Type[ProteinsComponentInput]:
        return ProteinsComponentInput

    @property
    def _output_parameter_type(self) -> Type[ProteinsComponentOutput]:
        return ProteinsComponentOutput
