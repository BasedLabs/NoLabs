import uuid
from typing import Optional

from nolabs.application.workflow import Component
from nolabs.application.workflow.models import ComponentDbModel


class WorkflowRepository:
    def fetch_component(self, component_id: uuid.UUID) -> Optional[Component]:
        model: ComponentDbModel = ComponentDbModel.objects.with_id(component_id)

        if not model:
            return None

        return model.get_component()

    def save_component(self, component: Component):
        model: ComponentDbModel = ComponentDbModel.objects.with_id(component.id)
        model.set_component(component=component)
        model.save()