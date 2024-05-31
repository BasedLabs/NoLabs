import pickle
import uuid
from typing import Dict, Any, List, Optional
from uuid import UUID

from mongoengine import Document, UUIDField, BinaryField, ReferenceField, CASCADE, DictField, StringField, IntField

from nolabs.refined.domain.models.common import Experiment, Job
from nolabs.workflow.workflow_schema import WorkflowSchemaModel


class WorkflowSchemaDbModel(Document):
    id: UUID = UUIDField(primary_key=True)
    experiment: Experiment = ReferenceField(Experiment, unique=True, required=True, reverse_delete_rule=CASCADE)
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


class PythonComponentDbModel(Document):
    id: uuid.UUID = UUIDField(primary_key=True)
    workflow: WorkflowSchemaDbModel = ReferenceField(WorkflowSchemaDbModel, reverse_delete_rule=CASCADE)

    input_parameter_dict: Dict[str, Any] = DictField(default=dict)
    output_parameter_dict: Dict[str, Any] = DictField(default=dict)

    jobs: List[Job] = ReferenceField(Job)

    def __init__(self, id: uuid.UUID, workflow: WorkflowSchemaDbModel, input_parameter_dict: Dict[str, Any],
                 output_parameter_dict: Dict[str, Any], jobs: List[Job], *args,
                 **values):
        self.id = id
        self.workflow = workflow
        self.input_parameter_dict = input_parameter_dict
        self.output_parameter_dict = output_parameter_dict
        self.jobs = jobs
        super().__init__(id=id,
                         workflow=workflow,
                         input_parameter_dict=input_parameter_dict,
                         output_parameter_dict=output_parameter_dict,
                         jobs=jobs,
                         *args, **values)
