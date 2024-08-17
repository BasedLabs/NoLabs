import logging
from typing import Optional
from uuid import UUID

from nolabs.application.use_cases.workflow.api_models import (GetComponentStateResponse, GetComponentStateRequest,
                                                              AllWorkflowDefinitionsResponse, ResetWorkflowRequest,
                                                              StartWorkflowComponentRequest)
from nolabs.application.workflow.api import Workflow, WorkflowDefinition
from nolabs.exceptions import NoLabsException, ErrorCodes


class DeleteWorkflowDefinitionFeature:
    async def handle(self, id: UUID):
        workflow = Workflow.get(id)
        workflow.delete()


class AllWorkflowDefinitionsFeature:
    async def handle(self, experiment_id: UUID) -> AllWorkflowDefinitionsResponse:
        ids = Workflow.all_workflow_ids(experiment_id=experiment_id)
        return AllWorkflowDefinitionsResponse(
            ids=ids
        )


class CreateWorkflowDefinitionFeature:
    async def handle(self, experiment_id: UUID) -> WorkflowDefinition:
        workflow = Workflow.create(experiment_id=experiment_id)
        return workflow.definition


class GetWorkflowDefinitionFeature:
    async def handle(self, id: UUID) -> Optional[WorkflowDefinition]:
        workflow = Workflow.get(id)

        if not workflow:
            raise NoLabsException(ErrorCodes.workflow_not_found)

        return workflow.definition


class UpdateWorkflowDefinitionFeature:
    async def handle(self, definition: WorkflowDefinition) -> WorkflowDefinition:
        workflow = Workflow.get(definition.workflow_id)
        workflow.update(definition=definition)
        return workflow.definition


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
