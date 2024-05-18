from abc import abstractmethod
from typing import TypeVar, Generic

from pydantic import BaseModel

TInput = TypeVar('TInput', bound=BaseModel)
TOutput = TypeVar('TOutput', bound=BaseModel)


class ComponentState(BaseModel):
    running: bool = False


class Component(Generic[TInput, TOutput]):
    async def state(self) -> ComponentState:
        return ComponentState(
            running=False
        )

    async def stop(self):
        pass

    @abstractmethod
    async def start(self, parameter: TInput) -> TOutput:
        ...

    @property
    def title(self) -> str:
        return 'Component'

    @property
    @abstractmethod
    def id(self) -> str:
        ...


