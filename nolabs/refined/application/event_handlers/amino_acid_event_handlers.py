__all__ = [
    'AminoAcidCreatedEventHandler'
]

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import AminoAcidCreatedEvent, AminoAcid
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class AminoAcidCreatedEventHandler(DomainEventHandler[AminoAcidCreatedEvent]):
    def handle(self, event: AminoAcidCreatedEvent):
        if AminoAcid.objects(name=event.amino_acid.name):
            raise NoLabsException('Cannot create amino acid with same name', error_code=ErrorCodes.duplicate_amino_acid)
