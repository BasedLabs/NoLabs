import uuid
from datetime import datetime
from typing import Any, Dict, List, TYPE_CHECKING, Optional, Union

from mongoengine import ReferenceField, UUIDField, Document, DictField, CASCADE, ListField, DateTimeField, StringField, \
    EmbeddedDocumentListField, IntField, EmbeddedDocument

if TYPE_CHECKING:
    from nolabs.application.workflow.component import Component


class PropertyErrorData(EmbeddedDocument):
    loc: List[str] = ListField(StringField())
    msg: str = StringField()

    @classmethod
    def create(cls, loc: List[str], msg: str) -> 'PropertyErrorData':
        return PropertyErrorData(
            loc=loc,
            msg=msg
        )


class WorkflowData(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    schema: Dict[str, Any] = DictField()

    last_executed_at: datetime = DateTimeField()
    flow_run_id: uuid.UUID = UUIDField()

    meta = {'collection': 'workflows'}

    @staticmethod
    def create(id: uuid.UUID,
               schema: Dict[str, Any]) -> 'WorkflowData':
        return WorkflowData(
            id=id,
            schema=schema
        )


class ComponentData(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    workflow: WorkflowData = ReferenceField(WorkflowData, reverse_delete_rule=CASCADE)

    input_errors: List[PropertyErrorData] = EmbeddedDocumentListField(PropertyErrorData, default=list)
    output_errors: List[PropertyErrorData] = EmbeddedDocumentListField(PropertyErrorData, default=list)

    executed_at: Optional[datetime] = DateTimeField()
    exception: Optional[str] = StringField()
    flow_run_id: Optional[uuid.UUID] = UUIDField()

    name: str = StringField()

    # region component fields

    input_schema: Dict[str, Any] = DictField()
    output_schema: Dict[str, Any] = DictField()
    input_value_dict: Dict[str, Any] = DictField()
    output_value_dict: Dict[str, Any] = DictField()
    previous_component_ids: List[uuid.UUID] = ListField(UUIDField())

    # endregion

    meta = {'collection': 'components'}

    @classmethod
    def create(cls,
               id: uuid.UUID,
               workflow: Union[WorkflowData, uuid.UUID]):
        return ComponentData(
            id=id,
            workflow=workflow
        )


class ExperimentWorkflowRelation(Document):
    id: str = StringField(primary_key=True)
    experiment = ReferenceField('Experiment', required=True, reverse_delete_rule=CASCADE)
    workflow: WorkflowData = ReferenceField(WorkflowData, required=True, reverse_delete_rule=CASCADE)

    @classmethod
    def create(cls, experiment_id: uuid.UUID, workflow_id: uuid.UUID) -> 'ExperimentWorkflowRelation':
        return ExperimentWorkflowRelation(id=f'{str(experiment_id)}|{str(workflow_id)}', experiment=experiment_id,
                                          workflow=workflow_id)


class JobRunData(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    task_run_id: uuid.UUID = UUIDField(required=True)
    executed_at: datetime = DateTimeField(default=datetime.utcnow, required=True)
    timeout: int = IntField(required=True)
    component: ComponentData = ReferenceField(ComponentData, required=True, reverse_delete_rule=CASCADE)
    component_id: uuid.UUID = UUIDField(required=True)

    @classmethod
    def create(cls,
               component_id: uuid.UUID,
               id: uuid.UUID,
               task_run_id: uuid.UUID,
               timeout: int,
               executed_at: datetime) -> 'JobRunData':
        return JobRunData(component=component_id, component_id=component_id, id=id, task_run_id=task_run_id,
                          timeout=timeout, executed_at=executed_at)
