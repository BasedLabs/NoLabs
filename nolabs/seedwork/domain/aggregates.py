from abc import abstractmethod
from dataclasses import field
from typing import Dict, Any

from pydantic import PrivateAttr, BaseModel

from nolabs.seedwork.domain.events import DomainEvent


class Aggregate:
    _events: list = PrivateAttr(default_factory=list)

    def register_event(self, event: DomainEvent):
        self._events.append(event)

    def collect_events(self):
        events = self._events
        self._events = []
        return events

    def clear_events(self):
        self._events = []