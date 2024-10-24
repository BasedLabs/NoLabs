from typing import Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends

from .api_models import AllWorkflowSchemasResponse, GetComponentStateRequest, ResetWorkflowRequest, \
    StartWorkflowComponentRequest, GetComponentStateResponse
from .di import WorkflowDependencies
from .use_cases import CreateWorkflowSchemaFeature, GetWorkflowSchemaFeature, UpdateWorkflowSchemaFeature, \
    StartWorkflowFeature, DeleteWorkflowSchemaFeature, AllWorkflowSchemasFeature, GetComponentStateFeature, \
    ResetWorkflowFeature, StartWorkflowComponentFeature
from nolabs.application.workflow.workflow_schema import WorkflowSchemaModel

router = APIRouter(
    prefix='/api/v1/workflow',
    tags=['Workflow'],
)


@router.post('/{experiment_id}', summary='Create workflow schema')
async def create_schema(
        feature: Annotated[
            CreateWorkflowSchemaFeature, Depends(WorkflowDependencies.create_workflow_schema)],
        experiment_id: UUID
) -> WorkflowSchemaModel:
    return await feature.handle(experiment_id=experiment_id)


@router.delete('/{workflow_id}', summary='Delete workflow schema')
async def create_schema(
        feature: Annotated[
            DeleteWorkflowSchemaFeature, Depends(WorkflowDependencies.delete_workflow_schema)],
        workflow_id: UUID
):
    return await feature.handle(workflow_id=workflow_id)


@router.get('/{workflow_id}', summary='Get workflow schema')
async def get_schema(
        feature: Annotated[
            GetWorkflowSchemaFeature, Depends(WorkflowDependencies.get_workflow_schema)],
        workflow_id: UUID
) -> Optional[WorkflowSchemaModel]:
    return await feature.handle(workflow_id=workflow_id)


@router.get('/all/{experiment_id}', summary='All workflow schemas')
async def get_schema(
        feature: Annotated[
            AllWorkflowSchemasFeature, Depends(WorkflowDependencies.all_workflow_schemas)],
        experiment_id: UUID
) -> AllWorkflowSchemasResponse:
    return await feature.handle(experiment_id=experiment_id)


@router.put('/', summary='Update workflow schema')
async def update_workflow_schema(
        feature: Annotated[
            UpdateWorkflowSchemaFeature, Depends(WorkflowDependencies.update_workflow_schema)],
        workflow_schema: WorkflowSchemaModel
) -> WorkflowSchemaModel:
    return await feature.handle(workflow_schema=workflow_schema)


@router.post('/{workflow_id}/start', summary='Start workflow')
async def start_workflow(
        feature: Annotated[
            StartWorkflowFeature, Depends(WorkflowDependencies.start_workflow)],
        workflow_id: UUID
):
    return await feature.handle(workflow_id=workflow_id)


@router.post('/{workflow_id}/start/{component_id}', summary='Start workflow component')
async def start_component(
        feature: Annotated[
            StartWorkflowComponentFeature, Depends(WorkflowDependencies.start_workflow_component)],
        workflow_id: UUID,
        component_id: UUID
):
    return await feature.handle(request=StartWorkflowComponentRequest(
        workflow_id=workflow_id,
        component_id=component_id
    ))


@router.post('/{workflow_id}/reset', summary='Reset workflow schema')
async def reset_workflow(
        feature: Annotated[
            ResetWorkflowFeature, Depends(WorkflowDependencies.reset_workflow)],
        request: ResetWorkflowRequest
):
    return await feature.handle(request=request)


@router.get('/component/{component_id}/state', summary='Get state')
async def get_component_state(
        feature: Annotated[
GetComponentStateFeature, Depends(WorkflowDependencies.get_component_state)
        ],
        component_id: UUID
) -> GetComponentStateResponse:
    return await feature.handle(request=GetComponentStateRequest(component_id=component_id))