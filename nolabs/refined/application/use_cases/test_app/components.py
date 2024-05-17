from typing import Optional

from pydantic import BaseModel

from nolabs.workflow.component import Component, TInput, TOutput


class Input(BaseModel):
    number: Optional[int] = 10


class Output(BaseModel):
    number: int


class PlusOneComponent(Component):
    @property
    def title(self) -> str:
        return 'One plus one component'

    @property
    def name(self) -> str:
        return 'OnePlusOneComponent'

    async def handle(self, parameter: Input) -> Output:
        return Output(
            number=parameter.number + 1  # type: ignore
        )


class PlusTwoComponent(Component):
    @property
    def title(self) -> str:
        return 'One plus one component'

    @property
    def name(self) -> str:
        return 'OnePlusOneComponent'

    async def handle(self, parameter: Input) -> Output:
        return Output(
            number=parameter.number + 2  # type: ignore
        )
