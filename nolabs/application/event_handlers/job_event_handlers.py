__all__ = [
    'JobStartedEventHandler',
    'JobFinishedEventHandler'
]

from nolabs.application.websockets.ws import send_websocket_message
from nolabs.domain.models.common import JobStartedEvent, JobFinishedEvent
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class JobStartedEventHandler(DomainEventHandler[JobStartedEvent]):
    async def handle(self, event: JobStartedEvent):
        await send_websocket_message(event_type='job_started', message={'job_id': event.job.iid.value})


class JobFinishedEventHandler(DomainEventHandler[JobFinishedEvent]):
    async def handle(self, event: JobFinishedEvent):
        await send_websocket_message(event_type='job_finished', message={'job_id': event.job.iid.value})
