from abc import abstractmethod, ABC
from typing import TypeVar, Generic

__all__ = ['EventHandler']


TEvent = TypeVar('TEvent')


class EventHandler(ABC, Generic[TEvent]):
    @abstractmethod
    def handle(self, event: TEvent):
        pass
