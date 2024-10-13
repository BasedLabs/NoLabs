from typing import Optional
from uuid import UUID

from fastapi import APIRouter

from workflow.application.schema import WorkflowSchema

from .api_models import (
    GetComponentRequest,
    GetComponentResponse,
    GetJobRequest,
    GetJobState,
    StartWorkflowComponentRequest,
)
from .use_cases import (
    CreateWorkflowSchemaFeature,
    GetComponentStateFeature,
    GetJobStateFeature,
    GetWorkflowSchemaFeature,
    StartWorkflowComponentFeature,
    StartWorkflowFeature,
    UpdateWorkflowSchemaFeature,
)

router = APIRouter(
    prefix="/api/v1/workflow",
    tags=["Workflow"],
)


@router.post("/{experiment_id}", summary="Create workflow schema")
async def create_workflow_schema(experiment_id: UUID,
) -> WorkflowSchema:
    return await CreateWorkflowSchemaFeature.handle(experiment_id=experiment_id)


@router.get("/{experiment_id}", summary="Get workflow schema")
async def get_schema() -> Optional[WorkflowSchema]:
    return await GetWorkflowSchemaFeature().handle(id=experiment_id)


@router.put("/", summary="Update workflow schema")
async def update_workflow_schema(schema: WorkflowSchema,
) -> WorkflowSchema:
    return await UpdateWorkflowSchemaFeature().handle(schema=schema)


@router.post("/{experiment_id}/start", summary="Start workflow")
async def start_workflow(experiment_id: UUID,
):
    return await StartWorkflowFeature().handle(id=experiment_id)


@router.post(
    "/{experiment_id}/start/{component_id}", summary="Start workflow component"
)
async def start_component(experiment_id: UUID, component_id: UUID):
    return await StartWorkflowComponentFeature().handle(
        request=StartWorkflowComponentRequest(
            experiment_id=experiment_id, component_id=component_id
        )
    )


@router.get("/job/{job_id}/state", summary="Get job state")
async def get_job_state(job_id: UUID) -> GetJobState:
    return await GetJobStateFeature().handle(request=GetJobRequest(job_id=job_id))


@router.get("/component/{component_id}/state", summary="Get state")
async def get_component_state(component_id: UUID) -> GetComponentResponse:
    return await GetComponentStateFeature().handle(request=GetComponentRequest(id=component_id))
