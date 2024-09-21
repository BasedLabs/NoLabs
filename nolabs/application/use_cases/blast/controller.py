__all__ = [
    "router",
]

from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.use_cases.blast.api_models import (
    GetJobStatusResponse, JobResponse, SetupJobRequest)
from nolabs.application.use_cases.blast.di import BlastDependencies
from nolabs.application.use_cases.blast.use_cases import (GetJobFeature,
                                                          RunJobFeature,
                                                          SetupJobFeature)

router = APIRouter(prefix="/api/v1/blast", tags=["Blast"])


@router.post("/jobs/run/{job_id}", summary="Start blast job")
async def start_job(
    feature: Annotated[RunJobFeature, Depends(BlastDependencies.run_job)], job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id)


@router.get("/jobs/{job_id}", summary="Get job")
async def get_job(
    job_id: UUID, feature: Annotated[GetJobFeature, Depends(BlastDependencies.get_job)]
) -> JobResponse:
    return await feature.handle(job_id=job_id)

@router.post("/jobs", summary="Setup job")
async def setup_job(
    request: SetupJobRequest,
    feature: Annotated[SetupJobFeature, Depends(BlastDependencies.setup_job)],
) -> JobResponse:
    return await feature.handle(request=request)
