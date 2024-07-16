from typing import Generic, List, get_args

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.seedwork.domain.event_handlers import DomainEventHandler, TDomainEvent


class EventDispatcher:
    event_handlers: List[DomainEventHandler[TDomainEvent]] = []

    @classmethod
    def add_event_handler(cls, event_handler: DomainEventHandler[TDomainEvent]):
        if not [eh for eh in cls.event_handlers if type(eh) == type(event_handler)]:
            cls.event_handlers.append(event_handler)

    @classmethod
    async def raise_event(cls, event: Generic[TDomainEvent]):
        event_handler: DomainEventHandler[TDomainEvent]
        for event_handler in [eh for eh in cls.event_handlers if get_args(eh.__class__.__orig_bases__[0])[0] == type(event)]:
            if not event_handler:
                raise NoLabsException(ErrorCodes.no_domain_event_handler)
            await event_handler.handle(event)
