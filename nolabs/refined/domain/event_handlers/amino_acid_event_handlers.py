from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.events.amino_acid_events import AminoAcidCreatedEvent
from nolabs.refined.domain.event_handlers import EventHandler


__all__ = ['AminoAcidCreatedEventHandler']

from nolabs.refined.domain.repository import Repository


class AminoAcidCreatedEventHandler(EventHandler[AminoAcidCreatedEvent]):
    _repository: Repository

    def __init__(self, repository: Repository):
        self._repository = repository

    def handle(self, event: AminoAcidCreatedEvent):
        amino_acid = self._repository.get_amino_acid(event.amino_acid_id)
        if self._repository.count_amino_acids(amino_acid.name) > 1:
            raise NoLabsException('Cannot add amino acid that already exists in system', ErrorCodes.duplicate_amino_acid)

