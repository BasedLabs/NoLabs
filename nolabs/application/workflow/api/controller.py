from typing import Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.workflow.api.schema import WorkflowSchema

from .api_models import (
    GetComponentRequest,
    GetComponentResponse,
    GetJobRequest,
    GetJobState,
    StartWorkflowComponentRequest,
)
from .di import WorkflowDependencies
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
async def create_workflow_schema(
    feature: Annotated[
        CreateWorkflowSchemaFeature,
        Depends(WorkflowDependencies.create_workflow_schema),
    ],
    experiment_id: UUID,
) -> WorkflowSchema:
    return await feature.handle(experiment_id=experiment_id)


@router.get("/{experiment_id}", summary="Get workflow schema")
async def get_schema(
    feature: Annotated[
        GetWorkflowSchemaFeature, Depends(WorkflowDependencies.get_workflow_schema)
    ],
    experiment_id: UUID,
) -> Optional[WorkflowSchema]:
    return await feature.handle(id=experiment_id)


@router.put("/", summary="Update workflow schema")
async def update_workflow_schema(
    feature: Annotated[
        UpdateWorkflowSchemaFeature,
        Depends(WorkflowDependencies.update_workflow_schema),
    ],
    schema: WorkflowSchema,
) -> WorkflowSchema:
    return await feature.handle(schema=schema)


@router.post("/{experiment_id}/start", summary="Start workflow")
async def start_workflow(
    feature: Annotated[
        StartWorkflowFeature, Depends(WorkflowDependencies.start_workflow)
    ],
    experiment_id: UUID,
):
    return await feature.handle(id=experiment_id)


@router.post(
    "/{experiment_id}/start/{component_id}", summary="Start workflow component"
)
async def start_component(
    feature: Annotated[
        StartWorkflowComponentFeature,
        Depends(WorkflowDependencies.start_workflow_component),
    ],
    experiment_id: UUID,
    component_id: UUID,
):
    return await feature.handle(
        request=StartWorkflowComponentRequest(
            experiment_id=experiment_id, component_id=component_id
        )
    )


@router.get("/job/{job_id}/state", summary="Get job state")
async def get_job_state(
    feature: Annotated[GetJobStateFeature, Depends(WorkflowDependencies.get_job_state)],
    job_id: UUID,
) -> GetJobState:
    return await feature.handle(request=GetJobRequest(job_id=job_id))


@router.get("/component/{component_id}/state", summary="Get state")
async def get_component_state(
    feature: Annotated[
        GetComponentStateFeature, Depends(WorkflowDependencies.get_component_state)
    ],
    component_id: UUID,
) -> GetComponentResponse:
    return await feature.handle(request=GetComponentRequest(id=component_id))
