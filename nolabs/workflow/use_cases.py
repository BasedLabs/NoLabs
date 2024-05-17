import uuid
from typing import Optional
from uuid import UUID

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Experiment
from nolabs.workflow.models import WorkflowSchema
from nolabs.workflow.mongoengine_models import WorkflowSchemaModel


class GetWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> Optional[WorkflowSchema]:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        workflow: WorkflowSchemaModel = WorkflowSchemaModel.objects(experiment=experiment)

        if not workflow:
            return None

        return workflow.get_workflow_value()


class SetWorkflowSchemaFeature:
    async def handle(self, workflow_schema: WorkflowSchema):
        experiment: Experiment = Experiment.objects.with_id(workflow_schema.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        workflow: WorkflowSchemaModel = WorkflowSchemaModel.objects(experiment=experiment)

        

        if not workflow:
            workflow = WorkflowSchemaModel.create(
                id=uuid.uuid4(),
                experiment=experiment,
                value=workflow_schema
            )
        else:
            workflow.set_workflow_value(value=workflow_schema)

        workflow.save(cascade=True)




