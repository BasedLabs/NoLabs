import uuid
from typing import Optional

from bson import ObjectId

from nolabs.application.workflow.component import Component
from nolabs.application.workflow.models import ComponentDbModel, WorkflowDbModel
from nolabs.domain.models.common import Experiment


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

    def fetch_experiment(self, component_id) -> Experiment:
        component_id = ObjectId(str(component_id))

        workflow: WorkflowDbModel = WorkflowDbModel.objects(components__in=[component_id]).first()

        return workflow.experiment