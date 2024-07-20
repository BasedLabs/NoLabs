from nolabs.application.event_handlers import ProteinCreatedEventHandler
from nolabs.application.event_handlers.job_event_handlers import JobStartedEventHandler, JobFinishedEventHandler
from nolabs.domain.event_dispatcher import EventDispatcher


class EventHandlersDependencies:
    @staticmethod
    def inject():
        event_handlers = [
            ProteinCreatedEventHandler(),
            JobStartedEventHandler(),
            JobFinishedEventHandler() # asd
        ]

        types = [type(eh) for eh in event_handlers]
        for event_handler in event_handlers:
            if types.count(type(event_handler)) > 1:
                raise RuntimeError(f'Event handler {type(event_handler)} was registered multiple times')

            EventDispatcher.add_event_handler(event_handler)