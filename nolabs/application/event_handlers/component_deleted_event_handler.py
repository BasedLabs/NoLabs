from nolabs.application.diffdock import DiffDockComponent

from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.workflow.core.component import ComponentDeletedEvent
from nolabs.domain.models.common import ComponentData
from nolabs.domain.models.diffdock import DiffDockBindingJob
from nolabs.seedwork.domain.event_handlers import DomainEventHandler


class ComponentDeletedEventHandler(DomainEventHandler[ComponentDeletedEvent]):
    async def handle(self, event: ComponentDeletedEvent):
        if event.component_data.name == DiffDockComponent.name:
            await self._diffdock_component_deleted(component_data=event.component_data)

    async def _diffdock_component_deleted(self, component_data: ComponentData):
        for job in DiffDockBindingJob.objects(
            component=component_data.id, celery_task_id__ne=None
        ).only("celery_task_id"):
            celery = get_celery_app()
            await celery.cancel_task(task_id=job.celery_task_id)
