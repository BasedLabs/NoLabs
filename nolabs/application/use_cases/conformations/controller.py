from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.use_cases.conformations.api_models import (
    GetJobStatusResponse, JobResponse, SetupJobRequest)
from nolabs.application.use_cases.conformations.di import \
    ConformationsDependencies
from nolabs.application.use_cases.conformations.use_cases import (
    GetJobFeature, GetJobStatusFeature, RunJobFeature, SetupJobFeature)

router = APIRouter(prefix="/api/v1/conformations", tags=["Conformations"])


@router.post("/jobs/run/{job_id}")
async def run_job(
    feature: Annotated[RunJobFeature, Depends(ConformationsDependencies.run_job)],
    job_id: UUID,
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get("/jobs/{job_id}")
async def get_job(
    job_id: UUID,
    feature: Annotated[GetJobFeature, Depends(ConformationsDependencies.get_job)],
) -> JobResponse:
    return await feature.handle(job_id)


@router.post("/jobs")
async def setup_job(
    request: SetupJobRequest,
    feature: Annotated[SetupJobFeature, Depends(ConformationsDependencies.setup_job)],
) -> JobResponse:
    return await feature.handle(request=request)


@router.get("/jobs/{job_id}/status", summary="Get job execution status")
async def get_job_status(
    job_id: UUID,
    feature: Annotated[
        GetJobStatusFeature, Depends(ConformationsDependencies.get_job_status)
    ],
) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)
