from typing import Generic

from dependency_injector.wiring import inject, Provide

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.event_handlers import EventHandlersContainer
from nolabs.seedwork.domain.event_handlers import DomainEventHandler, TDomainEvent


class EventDispatcher:
    @staticmethod
    @inject
    def raise_event(event: Generic[TDomainEvent],
                    container: EventHandlersContainer = Provide[EventHandlersContainer]):
        event_handler: DomainEventHandler[TDomainEvent]
        while ((event_handler := next(container.traverse(types=[DomainEventHandler[TDomainEvent]]), None))
               is not None):
            if not event_handler:
                raise NoLabsException(f'No event handler was registered for domain event {str(event)}',
                                      ErrorCodes.no_domain_event_handler)
            event_handler.handle(event)
