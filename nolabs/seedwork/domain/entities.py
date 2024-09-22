__all__ = ["Entity"]

from typing import List

from nolabs.seedwork.domain.events import DomainEvent


class Entity:
    _events: List[DomainEvent] = []

    def register_event(self, event: DomainEvent):
        self._events.append(event)

    def collect_events(self) -> List[DomainEvent]:
        events = self._events
        self._events = []
        return events

    def clear_events(self):
        self._events = []
