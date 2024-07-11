__all__ = [
    'LigandCreatedEventHandler'
]

from nolabs.domain.models.common import LigandCreatedEvent
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class LigandCreatedEventHandler(DomainEventHandler[LigandCreatedEvent]):
    def handle(self, event: LigandCreatedEvent):
        pass
