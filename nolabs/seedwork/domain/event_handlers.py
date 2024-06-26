from abc import abstractmethod
from typing import TypeVar, Generic

from nolabs.seedwork.domain.events import DomainEvent

TDomainEvent = TypeVar('TDomainEvent', bound=DomainEvent)


class DomainEventHandler(Generic[TDomainEvent]):
    @abstractmethod
    def handle(self, event: TDomainEvent):
        pass
