__all__ = [
    "router",
]

from uuid import UUID

from fastapi import APIRouter

from nolabs.application.folding.api_models import JobResponse, SetupJobRequest
from nolabs.application.folding.use_cases import (
    GetJobFeature,
    RunJobFeature,
    SetupJobFeature,
)

router = APIRouter(
    prefix="/api/v1/folding",
    tags=["Folding"],
)


@router.post("/jobs/run/{job_id}", summary="Start folding job")
async def start_job(job_id: UUID) -> JobResponse:
    return await RunJobFeature().handle(job_id)


@router.get("/jobs/{job_id}", summary="Get job")
async def get_job(job_id: UUID) -> JobResponse:
    return await GetJobFeature().handle(job_id=job_id)


@router.post("/jobs", summary="Setup job")
async def setup_job(request: SetupJobRequest) -> JobResponse:
    return await SetupJobFeature().handle(request=request)
