import uuid

from mongoengine import Document, ReferenceField, CASCADE, StringField

from nolabs.application.workflow import WorkflowData
from nolabs.domain.models.common import Experiment


class ExperimentWorkflowRelation(Document):
    id: str = StringField(primary_key=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    workflow: WorkflowData = ReferenceField(WorkflowData, required=True, reverse_delete_rule=CASCADE)

    @classmethod
    def create(cls, experiment_id: uuid.UUID, workflow_id: uuid.UUID) -> 'ExperimentWorkflowRelation':
        return ExperimentWorkflowRelation(id=f'{str(experiment_id)}|{workflow_id: uuid.UUID}', experiment=experiment_id, workflow=workflow_id)
