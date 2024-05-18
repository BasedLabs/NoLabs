import uuid
from uuid import UUID
import pickle

from mongoengine import Document, UUIDField, BinaryField, ReferenceField, CASCADE

from nolabs.refined.domain.models.common import Experiment
from nolabs.workflow.workflow_schema import WorkflowSchemaModel


class WorkflowSchemaDbModel(Document):
    id: UUID = UUIDField(primary_key=True, default=uuid.uuid4)
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
