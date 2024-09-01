from mongoengine import Document, ReferenceField, CASCADE, StringField

from nolabs.application.workflow.data import WorkflowState
from nolabs.domain.models.common import Experiment


class ExperimentWorkflowRelation(Document):
    id: str = StringField(primary_key=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    workflow: WorkflowState = ReferenceField(WorkflowState, required=True, reverse_delete_rule=CASCADE)

    @classmethod
    def create(cls, experiment: Experiment, workflow: WorkflowState) -> 'ExperimentWorkflowRelation':
        return ExperimentWorkflowRelation(id=f'{str(experiment.id)}|{workflow.id}', experiment=experiment, workflow=workflow)
