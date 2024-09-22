__all__ = [
    "router",
]

from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.use_cases.binding_pockets.api_models import (
    GetJobStatusResponse, JobResponse, SetupJobRequest)
from nolabs.application.use_cases.binding_pockets.di import \
    BindingPocketsDependencies
from nolabs.application.use_cases.binding_pockets.use_cases import (
    GetJobStatusFeature, RunJobFeature, SetupJobFeature)

router = APIRouter(
    prefix="/api/v1/binding-pockets",
    tags=["Binding pockets"],
)


@router.post("/jobs/run/{job_id}", summary="Start binding pockets prediction job")
async def run_job(
    feature: Annotated[RunJobFeature, Depends(BindingPocketsDependencies.run_job)],
    job_id: UUID,
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get("/jobs/{job_id}", summary="Get job")
async def get_job(
    job_id: UUID,
    feature: Annotated[RunJobFeature, Depends(BindingPocketsDependencies.get_job)],
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.post("/jobs", summary="Setup job")
async def setup_job(
    request: SetupJobRequest,
    feature: Annotated[SetupJobFeature, Depends(BindingPocketsDependencies.setup_job)],
) -> JobResponse:
    return await feature.handle(request=request)


@router.get("/jobs/{job_id}/status", summary="Job status")
async def get_job_status(
    job_id: UUID,
    feature: Annotated[
        GetJobStatusFeature, Depends(BindingPocketsDependencies.get_job_status)
    ],
) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)
