__all__ = [
    'ProteinCreatedEventHandler'
]

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models import ProteinCreatedEvent, Protein
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class ProteinCreatedEventHandler(DomainEventHandler[ProteinCreatedEvent]):
    def handle(self, event: ProteinCreatedEvent):
        if Protein.objects(name=event.protein.name.value):
            raise NoLabsException(error_code=ErrorCodes.duplicate_protein)
