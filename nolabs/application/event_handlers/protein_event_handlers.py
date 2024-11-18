__all__ = ["ProteinCreatedEventHandler"]

from nolabs.domain.models.common import ProteinCreatedEvent
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class ProteinCreatedEventHandler(DomainEventHandler[ProteinCreatedEvent]):
    def handle(self, event: ProteinCreatedEvent):
        pass
