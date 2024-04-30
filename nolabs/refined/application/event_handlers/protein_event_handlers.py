__all__ = [
    'ProteinCreatedEventHandler'
]

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models import ProteinCreatedEvent, Protein
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class ProteinCreatedEventHandler(DomainEventHandler[ProteinCreatedEvent]):
    def handle(self, event: ProteinCreatedEvent):
        if Protein.objects(name=event.protein.name):
            raise NoLabsException('Cannot create protein with same name', error_code=ErrorCodes.duplicate_protein)
