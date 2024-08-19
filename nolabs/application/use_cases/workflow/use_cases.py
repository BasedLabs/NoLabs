import logging
from typing import Optional
from uuid import UUID

from nolabs.application.use_cases.workflow.api_models import (GetComponentStateResponse, GetComponentStateRequest,
                                                              AllWorkflowSchemasResponse, ResetWorkflowRequest,
                                                              StartWorkflowComponentRequest)
from nolabs.application.workflow import Workflow, WorkflowSchema
from nolabs.exceptions import NoLabsException, ErrorCodes


class DeleteWorkflowSchemaFeature:
    async def handle(self, id: UUID):
        workflow = Workflow.get(id)
        workflow.delete()


class AllWorkflowSchemasFeature:
    async def handle(self, experiment_id: UUID) -> AllWorkflowSchemasResponse:
        ids = Workflow.all_workflow_ids(experiment_id=experiment_id)
        return AllWorkflowSchemasResponse(
            ids=ids
        )


class CreateWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> WorkflowSchema:
        workflow = Workflow.create(experiment_id=experiment_id)
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
        await workflow.start()


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
