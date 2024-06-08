import pickle
import uuid
from typing import Dict, Any, List, Optional
from uuid import UUID

from mongoengine import Document, UUIDField, BinaryField, ReferenceField, CASCADE, DictField, StringField, IntField, \
    ListField

from nolabs.refined.domain.models.common import Experiment, Job
from nolabs.workflow.workflow_schema import WorkflowSchemaModel


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


class ComponentDbModel(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    workflow: WorkflowSchemaDbModel = ReferenceField(WorkflowSchemaDbModel, reverse_delete_rule=CASCADE)

    last_exception: Optional[str] = StringField(required=False)

    input_parameter_dict: Dict[str, Any] = DictField(default=dict)
    output_parameter_dict: Dict[str, Any] = DictField(default=dict)

    jobs: List[Job] = ListField(ReferenceField(Job), default=[])

    @classmethod
    def create(cls, id: uuid.UUID, workflow: WorkflowSchemaDbModel,
               last_exception: Optional[str], input_parameter_dict: Dict[str, Any],
                 output_parameter_dict: Dict[str, Any], jobs: List[Job]):
        return ComponentDbModel(
            id=id,
            last_exception=last_exception,
            workflow=workflow,
            input_parameter_dict=input_parameter_dict,
            output_parameter_dict=output_parameter_dict,
            jobs=jobs
        )
