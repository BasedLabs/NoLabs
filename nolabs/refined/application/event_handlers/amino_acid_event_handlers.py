__all__ = [
    'AminoAcidCreatedEvent',
    'AminoAcidCreatedEventHandler'
]

from pydantic.dataclasses import dataclass

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models import AminoAcid
from nolabs.refined.domain.repository import Repository
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


@dataclass
class AminoAcidCreatedEvent:
    amino_acid: AminoAcid


class AminoAcidCreatedEventHandler(DomainEventHandler[AminoAcidCreatedEvent]):
    _repository: Repository

    def __init__(self, repository: Repository):
        self._repository = repository

    def handle(self, event: AminoAcidCreatedEvent):
        if self._repository.amino_acids(name=event.amino_acid.name):
            raise NoLabsException('Cannot create amino acid with same name', error_code=ErrorCodes.duplicate_amino_acid)
