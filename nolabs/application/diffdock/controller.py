__all__ = [
    "router",
]

from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.diffdock.api_models import JobResponse, SetupJobRequest
from nolabs.application.diffdock.di import DiffDockDependencies
from nolabs.application.diffdock.use_cases import (
    GetJobFeature,
    RunJobFeature,
    SetupJobFeature,
)

router = APIRouter(prefix="/api/v1/diffdock", tags=["Diffdock"])


@router.post("/jobs/run/{job_id}", summary="Start diffdock job")
async def start_job(
    feature: Annotated[RunJobFeature, Depends(DiffDockDependencies.run_job)],
    job_id: UUID,
) -> JobResponse:
    return await feature.handle(job_id)


@router.get("/jobs/{job_id}", summary="Get job")
async def get_job(
    job_id: UUID,
    feature: Annotated[GetJobFeature, Depends(DiffDockDependencies.get_job)],
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.post("/jobs", summary="Setup job")
async def setup_job(
    request: SetupJobRequest,
    feature: Annotated[SetupJobFeature, Depends(DiffDockDependencies.setup_job)],
) -> JobResponse:
    return await feature.handle(request=request)
