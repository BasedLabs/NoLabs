import uuid
from datetime import datetime
from typing import Any, Dict, List, TYPE_CHECKING, Optional

from mongoengine import ReferenceField, UUIDField, Document, DictField, CASCADE, ListField, DateTimeField, StringField, \
    EmbeddedDocumentListField, IntField, EmbeddedDocument

from nolabs.application.workflow.schema import WorkflowSchema

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


class JobRunData(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    task_run_id: uuid.UUID = UUIDField(required=True)
    executed_at: datetime = DateTimeField(default=datetime.utcnow, required=True)
    timeout: int = IntField(required=True)

    state: str = StringField()
    state_message = StringField()

    component: 'ComponentData' = ReferenceField('ComponentData')

    @classmethod
    def create(cls,
               component_id: uuid.UUID,
               id: uuid.UUID,
               task_run_id: uuid.UUID,
               timeout: int,
               state: str,
               executed_at: datetime) -> 'JobRunData':
        return JobRunData(component=component_id, id=id, task_run_id=task_run_id, timeout=timeout, state=state, executed_at=executed_at)


class WorkflowData(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    schema: Dict[str, Any] = DictField()

    last_executed_at: datetime = DateTimeField()
    flow_run_id: uuid.UUID = UUIDField()

    state: str = StringField()
    state_message: str = StringField()

    meta = {'collection': 'workflows'}

    @staticmethod
    def create(id: uuid.UUID,
               schema: WorkflowSchema) -> 'WorkflowData':
        return WorkflowData(
            id=id,
            schema=schema.dict()
        )

    def set_schema(self, schema: WorkflowSchema):
        self.schema = schema.dict()

    def get_schema(self) -> WorkflowSchema:
        return WorkflowSchema(**self.schema)


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

    state: str = StringField()
    state_message: str = StringField()

    # endregion

    meta = {'collection': 'components'}

    @classmethod
    def create(cls,
               id: uuid.UUID,
               workflow: WorkflowData,
               component: 'Component'):
        return ComponentData(
            id=id,
            workflow=workflow,
            input_schema=component.input_schema.dict(),
            output_schema=component.output_schema.dict(),
            name=component.name
        )
