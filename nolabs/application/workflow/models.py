import uuid
from datetime import datetime
from typing import Dict, Any, List
from uuid import UUID

from mongoengine import Document, UUIDField, ReferenceField, CASCADE, DictField, StringField, IntField, \
    ListField, EmbeddedDocument, EmbeddedDocumentListField, PULL, DateTimeField

from nolabs.application.workflow.component import Component
from nolabs.application.workflow.definition import WorkflowDefinition
from nolabs.domain.models.common import Experiment, Job


class WorkflowDbModel(Document):
    id: UUID = UUIDField(primary_key=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    components: List['ComponentDbModel'] = ListField(ReferenceField(Job, reverse_delete_rule=PULL), default=[])
    error: str = StringField()
    definition: Dict[str, Any] = DictField()

    @staticmethod
    def create(id: UUID,
               experiment: Experiment,
               components: List['ComponentDbModel'],
               definition: WorkflowDefinition) -> 'WorkflowDbModel':
        return WorkflowDbModel(
            id=id,
            experiment=experiment,
            components=components,
            definition=definition.dict()
        )

    def get_workflow_definition(self) -> 'WorkflowDefinition':
        return WorkflowDefinition(**self.definition)

    def set_workflow_definition(self, value: WorkflowDefinition):
        self.definition = value.dict()


class InputPropertyErrorDbModel(EmbeddedDocument):
    loc: List[str] = ListField(StringField())
    msg: str = StringField()

    @classmethod
    def create(cls, loc: List[str], msg: str) -> 'InputPropertyErrorDbModel':
        return InputPropertyErrorDbModel(
            loc=loc,
            msg=msg
        )


class ComponentDbModel(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    workflow: WorkflowDbModel = ReferenceField(WorkflowDbModel, reverse_delete_rule=CASCADE)

    input_property_errors: List[InputPropertyErrorDbModel] = EmbeddedDocumentListField(InputPropertyErrorDbModel)

    jobs: List[Job] = ListField(ReferenceField(Job, reverse_delete_rule=PULL), default=[])
    last_jobs_count: int = IntField()

    last_executed_at: datetime = DateTimeField()

    component: Dict[str, Any] = DictField()

    @classmethod
    def create(cls,
               id: uuid.UUID,
               workflow: WorkflowDbModel,
               component: Component,
               jobs: List[Job]):
        return ComponentDbModel(
            id=id,
            workflow=workflow,
            component=component.dict(),
            jobs=jobs
        )

    def get_component(self) -> Component:
        return Component(**self.component)

    def set_component(self, component: Component):
        self.component = component.dict()

