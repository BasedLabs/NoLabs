import uuid
from datetime import datetime
from typing import Any, Dict, List, TYPE_CHECKING

from mongoengine import ReferenceField, UUIDField, Document, DictField, CASCADE, ListField, DateTimeField, StringField, \
    EmbeddedDocumentListField, IntField, EmbeddedDocument, EmbeddedDocumentField

from nolabs.application.workflow.schema import WorkflowSchema

if TYPE_CHECKING:
    from nolabs.application.workflow.component import Component


class InputPropertyErrorDbModel(EmbeddedDocument):
    loc: List[str] = ListField(StringField())
    msg: str = StringField()

    @classmethod
    def create(cls, loc: List[str], msg: str) -> 'InputPropertyErrorDbModel':
        return InputPropertyErrorDbModel(
            loc=loc,
            msg=msg
        )


class JobRunModel(EmbeddedDocument):
    id: uuid.UUID = UUIDField(primary_key=True)
    task_run_id: uuid.UUID = UUIDField(required=True)
    created_at: datetime = DateTimeField(default=datetime.utcnow, required=True)

    @classmethod
    def create(cls, id: uuid.UUID, task_run_id: uuid.UUID, created_at: datetime) -> 'JobRunModel':
        return JobRunModel(id=id, task_run_id=task_run_id, created_at=created_at)


class ComponentRunModel(EmbeddedDocument):
    task_run_id: uuid.UUID = UUIDField(primary_key=True, default=uuid.uuid4)
    created_at: datetime = DateTimeField(default=datetime.utcnow, required=True)
    jobs: List[JobRunModel] = EmbeddedDocumentListField(JobRunModel, default=[])

    @classmethod
    def create(cls, task_run_id: uuid.UUID, created_at: datetime) -> 'ComponentRunModel':
        return ComponentRunModel(task_run_id=task_run_id, created_at=created_at)


class WorkflowRunModel(EmbeddedDocument):
    flow_run_id: uuid.UUID = UUIDField(required=True)
    created_at: datetime = DateTimeField(default=datetime.utcnow, required=True)


class WorkflowState(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    schema: Dict[str, Any] = DictField()

    run: WorkflowRunModel = EmbeddedDocumentField(WorkflowRunModel)

    meta = {'collection': 'workflows'}

    @staticmethod
    def create(id: uuid.UUID,
               schema: WorkflowSchema) -> 'WorkflowState':
        return WorkflowState(
            id=id,
            schema=schema.dict()
        )

    def set_schema(self, schema: WorkflowSchema):
        self.schema = schema.dict()

    def get_schema(self) -> WorkflowSchema:
        return WorkflowSchema(**self.schema)


class ComponentState(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    workflow: WorkflowState = ReferenceField(WorkflowState, reverse_delete_rule=CASCADE)

    input_property_errors: List[InputPropertyErrorDbModel] = EmbeddedDocumentListField(InputPropertyErrorDbModel)

    last_jobs_count: int = IntField()

    last_executed_at: datetime = DateTimeField()

    name: str = StringField()

    # region component fields

    input_schema: Dict[str, Any] = DictField()
    output_schema: Dict[str, Any] = DictField()
    input_value_dict: Dict[str, Any] = DictField()
    output_value_dict: Dict[str, Any] = DictField()
    previous_component_ids: List[uuid.UUID] = ListField(UUIDField())

    runs: List[ComponentRunModel] = EmbeddedDocumentListField(ComponentRunModel)

    # endregion

    meta = {'collection': 'components'}

    @classmethod
    def create(cls,
               id: uuid.UUID,
               workflow: WorkflowState,
               component: 'Component'):
        return ComponentState(
            id=id,
            workflow=workflow,
            input_schema=component.input_schema.dict(),
            output_schema=component.output_schema.dict(),
            input_value_dict=component.input_value_dict,
            output_value_dict=component.output_value_dict,
            previous_component_ids=component.previous_component_ids,
            name=component.name,
            runs=[]
        )

    def set_component(self, component: 'Component'):
        self.input_schema = component.input_schema.dict()
        self.output_schema = component.output_schema.dict()
        self.input_value_dict = component.input_value_dict
        self.output_value_dict = component.output_value_dict
        self.previous_component_ids = component.previous_component_ids

    def latest_job_ids(self) -> List[uuid.UUID]:
        if not self.runs:
            return []

        run = self.runs[-1]
        return [j.id for j in run.jobs]