__all__ = [
    'ProteinCreatedEvent',
    'ProteinCreatedEventHandler'
]

from pydantic.dataclasses import dataclass

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models import Protein
from nolabs.refined.domain.repository import Repository
from nolabs.seedwork.domain.event_handlers import DomainEventHandler
from nolabs.seedwork.domain.events import DomainEvent


@dataclass
class ProteinCreatedEvent(DomainEvent):
    protein: Protein


class ProteinCreatedEventHandler(DomainEventHandler[ProteinCreatedEvent]):
    _repository: Repository

    def __init__(self, repository: Repository):
        self._repository = repository

    def handle(self, event: ProteinCreatedEvent):
        if self._repository.proteins(name=event.protein.name):
            raise NoLabsException('Cannot create protein with same name', error_code=ErrorCodes.duplicate_protein)
