import uuid
from typing import List, Optional, Type

from pydantic import BaseModel

from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler


class LigandsComponentInput(BaseModel):
    ligands: List[uuid.UUID]


class LigandsComponentOutput(BaseModel):
    ligands: List[uuid.UUID]


class LigandsComponent(Component[LigandsComponentInput, LigandsComponentOutput]):
    @property
    def output_parameter_type(self) -> Type[LigandsComponentOutput]:
        return LigandsComponentOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlowHandler"]:
        return LigandsFlow

    @property
    def input_parameter_type(self) -> Type[LigandsComponentInput]:
        return LigandsComponentInput

    name = "Ligands"
    description = "Ligands datasource"


class LigandsFlow(ComponentFlowHandler):
    async def on_finish(
        self, inp: LigandsComponentInput, job_ids: List[uuid.UUID]
    ) -> Optional[LigandsComponentOutput]:
        return LigandsComponentOutput(ligands=inp.ligands)
