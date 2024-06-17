__all__ = [
    'LigandCreatedEventHandler'
]

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import LigandCreatedEvent, Ligand
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class LigandCreatedEventHandler(DomainEventHandler[LigandCreatedEvent]):
    def handle(self, event: LigandCreatedEvent):
        pass
