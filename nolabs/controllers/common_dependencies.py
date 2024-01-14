from nolabs.features.events_queue import EventsQueue
from nolabs.infrastructure.settings import Settings


def settings_dependency() -> Settings:
    return Settings()


_events_queue = EventsQueue()


def events_queue_dependency() -> EventsQueue:
    return _events_queue
