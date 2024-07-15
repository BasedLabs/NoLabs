__all__ = [
    'JobStartedEventHandler',
    'JobFinishedEventHandler'
]

from nolabs.domain.models.common import JobStartedEvent, JobFinishedEvent
from nolabs.infrastructure.websocket_queue import websockets_queue
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class JobStartedEventHandler(DomainEventHandler[JobStartedEvent]):
    async def handle(self, event: JobStartedEvent):
        websockets_queue.write({'job_id': str(event.job.iid.value), 'event_type': 'job_started'})


class JobFinishedEventHandler(DomainEventHandler[JobFinishedEvent]):
    async def handle(self, event: JobFinishedEvent):
        websockets_queue.write({'job_id': str(event.job.iid.value), 'event_type': 'job_finished'})
