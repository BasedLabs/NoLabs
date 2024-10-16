import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.application.workflow.component import Component, JobValidationError


class LigandsComponentInput(BaseModel):
    ligands: List[uuid.UUID]


class LigandsComponentOutput(BaseModel):
    ligands: List[uuid.UUID]


class LigandsComponent(Component[LigandsComponentInput, LigandsComponentOutput]):
    name = 'Ligands'
    description = 'Ligands datasource'

    async def execute(self):
        self.output = LigandsComponentOutput(
            ligands=self.input.ligands
        )

    async def setup_jobs(self):
        pass

    async def jobs_setup_errors(self) -> List[JobValidationError]:
        return []

    @property
    def _input_parameter_type(self) -> Type[LigandsComponentInput]:
        return LigandsComponentInput

    @property
    def _output_parameter_type(self) -> Type[LigandsComponentOutput]:
        return LigandsComponentOutput
