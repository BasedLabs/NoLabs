from nolabs.application.diffdock import DiffDockComponent
from workflow.core.component import ComponentDeletedEvent
from domain.models.common import ComponentData
from nolabs.domain.models.diffdock import DiffDockBindingJob
from nolabs.infrastructure.cel import cel as celery
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class ComponentDeletedEventHandler(DomainEventHandler[ComponentDeletedEvent]):
    async def handle(self, event: ComponentDeletedEvent):
        if event.component_data.name == DiffDockComponent.name:
            await self._diffdock_component_deleted(component_data=event.component_data)

    async def _diffdock_component_deleted(self, component_data: ComponentData):
        for job in DiffDockBindingJob.objects(
            component=component_data.id, celery_task_id__ne=None
        ).only("celery_task_id"):
            await celery.cancel_task(task_id=job.celery_task_id)
