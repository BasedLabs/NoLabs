import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.application.workflow import Component, ComponentFlow
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
    def component_flow_type(self) -> Type['ComponentFlow']:
        return ProteinsTask


class ProteinsTask(ComponentFlow):
    component_timeout_seconds = 1.0

    async def post_execute(self, inp: ProteinsComponentInput, **kwargs):
        return ProteinsComponentOutput(proteins=inp.proteins)
