__all__ = [
    'ProteinCreatedEventHandler'
]

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import ProteinCreatedEvent, Protein
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class ProteinCreatedEventHandler(DomainEventHandler[ProteinCreatedEvent]):
    def handle(self, event: ProteinCreatedEvent):
        pass
