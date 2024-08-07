import pickle
import uuid
from datetime import datetime
from typing import Dict, Any, List
from uuid import UUID

from mongoengine import Document, UUIDField, BinaryField, ReferenceField, CASCADE, DictField, StringField, IntField, \
    ListField, EmbeddedDocument, EmbeddedDocumentListField, PULL, DateTimeField

from nolabs.domain.models.common import Experiment, Job
from nolabs.application.workflow.workflow_schema import WorkflowSchemaModel


class WorkflowSchemaDbModel(Document):
    id: UUID = UUIDField(primary_key=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    value: bytes = BinaryField(required=True)

    @staticmethod
    def create(id: UUID,
               experiment: Experiment,
               value: WorkflowSchemaModel) -> 'WorkflowSchemaDbModel':
        bytes_value = pickle.dumps(value)

        return WorkflowSchemaDbModel(
            id=id,
            experiment=experiment,
            value=bytes_value
        )

    def get_workflow_value(self) -> 'WorkflowSchemaModel':
        return pickle.loads(self.value)

    def set_workflow_value(self, value: WorkflowSchemaModel):
        self.value = pickle.dumps(value)


class JobErrorDbModel(EmbeddedDocument):
    job_id: UUID = UUIDField()
    msg: str = StringField()

    @classmethod
    def create(cls, job_id: UUID, msg: str) -> 'JobErrorDbModel':
        return JobErrorDbModel(
            job_id=job_id,
            msg=msg
        )


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
    workflow: WorkflowSchemaDbModel = ReferenceField(WorkflowSchemaDbModel, reverse_delete_rule=CASCADE)

    last_exceptions: List[str] = ListField(StringField(required=False))
    input_property_errors: List[InputPropertyErrorDbModel] = EmbeddedDocumentListField(InputPropertyErrorDbModel)

    input_parameter_dict: Dict[str, Any] = DictField(default=dict)
    output_parameter_dict: Dict[str, Any] = DictField(default=dict)

    jobs: List[Job] = ListField(ReferenceField(Job, reverse_delete_rule=PULL), default=[])
    last_jobs_count: int = IntField()

    last_executed_at: datetime = DateTimeField()

    @classmethod
    def create(cls,
               id: uuid.UUID,
               workflow: WorkflowSchemaDbModel,
               input_parameter_dict: Dict[str, Any],
               output_parameter_dict: Dict[str, Any],
               jobs: List[Job],
               input_property_errors: List[InputPropertyErrorDbModel],
               last_exceptions: List[str]):
        return ComponentDbModel(
            id=id,
            workflow=workflow,
            input_parameter_dict=input_parameter_dict,
            output_parameter_dict=output_parameter_dict,
            input_property_errors=input_property_errors,
            jobs=jobs,
            last_exceptions=last_exceptions
        )
