from abc import abstractmethod
from typing import Generic, TypeVar

from nolabs.seedwork.domain.events import DomainEvent

TDomainEvent = TypeVar("TDomainEvent", bound=DomainEvent)


class DomainEventHandler(Generic[TDomainEvent]):
    @abstractmethod
    def handle(self, event: TDomainEvent):
        pass
