import logging
import uuid
from typing import Optional
from uuid import UUID

from nolabs.application.use_cases.workflow.api_models import (GetComponentStateResponse, GetComponentStateRequest,
                                                              AllWorkflowSchemasResponse, ResetWorkflowRequest,
                                                              StartWorkflowComponentRequest)
from nolabs.application.workflow import Workflow, WorkflowSchema
from nolabs.domain.models.common import Experiment
from nolabs.exceptions import NoLabsException, ErrorCodes


class DeleteWorkflowSchemaFeature:
    async def handle(self, id: UUID):
        workflow = Workflow.get(id)
        workflow.delete()


class AllWorkflowSchemasFeature:
    async def handle(self, experiment_id: UUID) -> AllWorkflowSchemasResponse:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)
        return AllWorkflowSchemasResponse(
            ids=experiment.workflow_ids
        )


class CreateWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> WorkflowSchema:
        workflow = Workflow.create(uuid.uuid4())
        experiment = Experiment.objects.with_id(experiment_id)
        experiment.add_workflow(workflow.id)
        return workflow.schema


class GetWorkflowSchemaFeature:
    async def handle(self, id: UUID) -> Optional[WorkflowSchema]:
        workflow = Workflow.get(id)

        if not workflow:
            raise NoLabsException(ErrorCodes.workflow_not_found)

        return workflow.schema


class UpdateWorkflowSchemaFeature:
    async def handle(self, schema: WorkflowSchema) -> WorkflowSchema:
        workflow = Workflow.get(schema.workflow_id)
        workflow.update(schema=schema)
        return workflow.schema


class StartWorkflowFeature:
    async def handle(self, id: UUID):
        workflow = Workflow.get(id)
        experiment = Experiment.objects(workflow_ids=id).first()
        await workflow.start(extra_dag_parameters={
            'experiment_id': experiment.id
        })


class StartWorkflowComponentFeature:
    logger: logging.Logger

    def __init__(self, logger: logging.Logger = logging.getLogger('nolabs')):
        self.logger = logger

    async def handle(self, request: StartWorkflowComponentRequest):
        pass


class GetComponentStateFeature:
    async def handle(self, request: GetComponentStateRequest) -> GetComponentStateResponse:
        pass


class ResetWorkflowFeature:
    async def handle(self, request: ResetWorkflowRequest):
        pass
