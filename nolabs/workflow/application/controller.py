from typing import Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends

from .api_models import AllWorkflowSchemasResponse
from .di import WorkflowDependencies
from .use_cases import CreateWorkflowSchemaFeature, GetWorkflowSchemaFeature, SetWorkflowSchemaFeature, \
    StartWorkflowFeature, DeleteWorkflowSchemaFeature, AllWorkflowSchemasFeature
from ..workflow_schema import WorkflowSchemaModel

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


@router.put('/', summary='Set workflow schema')
async def set_workflow_schema(
        feature: Annotated[
            SetWorkflowSchemaFeature, Depends(WorkflowDependencies.set_workflow_schema)],
        workflow_schema: WorkflowSchemaModel
) -> WorkflowSchemaModel:
    return await feature.handle(workflow_schema=workflow_schema)


@router.post('/{experiment_id}/start', summary='Start workflow schema')
async def start_workflow(
        feature: Annotated[
            StartWorkflowFeature, Depends(WorkflowDependencies.start_workflow)],
        experiment_id: UUID
):
    return await feature.handle(experiment_id=experiment_id)
