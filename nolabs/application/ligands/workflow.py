import uuid
from typing import List, Optional, Type

from pydantic import BaseModel

from workflow.flows import ComponentFlow
from nolabs.domain.workflow.component import Component, TInput, TOutput


class LigandsComponentInput(BaseModel):
    ligands: List[uuid.UUID]


class LigandsComponentOutput(BaseModel):
    ligands: List[uuid.UUID]


class LigandsComponent(Component[LigandsComponentInput, LigandsComponentOutput]):
    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return LigandsComponentOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlow"]:
        return LigandsFlow

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return LigandsComponentInput

    name = "Ligands"
    description = "Ligands datasource"


class LigandsFlow(ComponentFlow):
    async def gather_jobs(
        self, inp: LigandsComponentInput, job_ids: List[uuid.UUID]
    ) -> Optional[TOutput]:
        return LigandsComponentOutput(ligands=inp.ligands)
