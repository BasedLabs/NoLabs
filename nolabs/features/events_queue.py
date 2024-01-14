from enum import Enum
from typing import Dict, Any


class EventsQueueMessageClass(Enum):
    conformations = 1


class EventsQueue:
    def __init__(self):
        self._events = {}

    def send_json(self, message_class: EventsQueueMessageClass, data: Dict[Any, Any]):
        self._events[message_class] = data

    def get_json(self, message_class: EventsQueueMessageClass) -> Dict[Any, Any] | None:
        return self._events.get(message_class)
