import uuid
from typing import List, Optional, Type

from pydantic import BaseModel

from nolabs.workflow.core.component import Component, TInput, TOutput
from nolabs.workflow.core.flow import ComponentFlowHandler


class ProteinsComponentInput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinsComponentOutput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinsComponent(Component[ProteinsComponentInput, ProteinsComponentOutput]):
    name = "Proteins"
    description = "Proteins datasource"

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return ProteinsComponentInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return ProteinsComponentOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlowHandler"]:
        return ProteinsFlow


class ProteinsFlow(ComponentFlowHandler):
    async def on_completion(
        self, inp: ProteinsComponentInput, job_ids: List[uuid.UUID]
    ) -> Optional[ProteinsComponentOutput]:
        return ProteinsComponentOutput(proteins=inp.proteins)
