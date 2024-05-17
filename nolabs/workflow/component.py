from abc import abstractmethod
from typing import TypeVar, Generic

from pydantic import BaseModel

TInput = TypeVar('TInput', bound=BaseModel)
TOutput = TypeVar('TOutput', bound=BaseModel)


class Component(Generic[TInput, TOutput]):
    @abstractmethod
    async def handle(self, parameter: TInput) -> TOutput:
        ...

    @property
    def title(self) -> str:
        return 'Component'

    @property
    @abstractmethod
    def id(self) -> str:
        ...


