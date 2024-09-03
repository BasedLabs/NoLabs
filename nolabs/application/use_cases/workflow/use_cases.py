import logging
import uuid
from typing import Optional, List
from uuid import UUID

from nolabs.application.use_cases.workflow.api_models import (GetComponentStateResponse, GetComponentStateRequest,
                                                              AllWorkflowSchemasResponse, ResetWorkflowRequest,
                                                              StartWorkflowComponentRequest)
from nolabs.application.use_cases.workflow.data import ExperimentWorkflowRelation
from nolabs.application.workflow import Workflow, WorkflowSchema, Component
from nolabs.application.workflow.workflow import get_component
from nolabs.domain.models.common import Experiment
from nolabs.exceptions import NoLabsException, ErrorCodes


class DeleteWorkflowSchemaFeature:
    async def handle(self, id: UUID):
        workflow = Workflow.get(id)
        workflow.delete()


class AllWorkflowSchemasFeature:
    async def handle(self, experiment_id: UUID) -> AllWorkflowSchemasResponse:
        relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(experiment=experiment_id).first()

        return AllWorkflowSchemasResponse(
            ids=[relation.workflow.id] if relation else []
        )


class CreateWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> WorkflowSchema:
        workflow = Workflow.create(uuid.uuid4())
        experiment = Experiment.objects.with_id(experiment_id)
        ExperimentWorkflowRelation.create(
            experiment=experiment,
            workflow=workflow.state
        ).save()
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
        relation: ExperimentWorkflowRelation = ExperimentWorkflowRelation.objects(workflow=workflow.id).first()

        if not relation:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        await workflow.start(extra_dag_parameters={
            'experiment_id': relation.experiment.id
        })


class StartWorkflowComponentFeature:
    logger: logging.Logger

    def __init__(self, logger: logging.Logger = logging.getLogger('nolabs')):
        self.logger = logger

    async def handle(self, request: StartWorkflowComponentRequest):
        pass


class GetComponentStateFeature:
    async def handle(self, request: GetComponentStateRequest) -> GetComponentStateResponse:
        component = get_component(request.component_id)

        return GetComponentStateResponse(
            input_dict=component.input_value_dict,
            output_dict=component.output_value_dict,
            job_ids=component.,
            input_property_errors=[],
            last_exceptions=[],
            jobs_errors=[]
        )


class ResetWorkflowFeature:
    async def handle(self, request: ResetWorkflowRequest):
        pass
