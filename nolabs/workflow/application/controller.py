from typing import Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends
from .di import WorkflowDependencies
from .use_cases import CreateWorkflowSchemaFeature, GetWorkflowSchemaFeature, SetWorkflowSchemaFeature, \
    StartWorkflowFeature, StopWorkflowFeature, GetComponentJobIdsFeature
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


@router.get('/{experiment_id}', summary='Get workflow schema')
async def get_schema(
        feature: Annotated[
            GetWorkflowSchemaFeature, Depends(WorkflowDependencies.get_workflow_schema)],
        experiment_id: UUID
) -> Optional[WorkflowSchemaModel]:
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


@router.post('/{experiment_id}/stop', summary='Stop workflow schema')
async def stop_workflow(
        feature: Annotated[
            StopWorkflowFeature, Depends(WorkflowDependencies.stop_workflow)],
        experiment_id: UUID
):
    return await feature.handle(experiment_id=experiment_id)


@router.post('/components/job-ids', summary='Get component job ids')
async def get_component_job_ids(
        feature: Annotated[
            GetComponentJobIdsFeature, Depends(WorkflowDependencies.get_component_job_ids)],
        experiment_id: UUID
):
    return await feature.handle()