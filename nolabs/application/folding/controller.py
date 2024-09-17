__all__ = [
    "router",
]

from typing import Annotated
from uuid import UUID

from application.folding.api_models import (GetJobStatusResponse, JobResponse,
                                            SetupJobRequest)
from application.folding.di import FoldingDependencies
from application.folding.use_cases import (GetJobFeature, GetJobStatusFeature,
                                           RunJobFeature, SetupJobFeature)
from fastapi import APIRouter, Depends

router = APIRouter(
    prefix="/api/v1/folding",
    tags=["Folding"],
)


@router.post("/jobs/run/{job_id}", summary="Start folding job")
async def start_job(
    feature: Annotated[RunJobFeature, Depends(FoldingDependencies.run_job)],
    job_id: UUID,
) -> JobResponse:
    return await feature.handle(job_id)


@router.get("/jobs/{job_id}", summary="Get job")
async def get_job(
    job_id: UUID,
    feature: Annotated[GetJobFeature, Depends(FoldingDependencies.get_job)],
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get("/jobs/{job_id}/status", summary="Get job execution status")
async def get_job_status(
    job_id: UUID,
    feature: Annotated[
        GetJobStatusFeature, Depends(FoldingDependencies.get_job_status)
    ],
) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)


@router.post("/jobs", summary="Setup job")
async def setup_job(
    request: SetupJobRequest,
    feature: Annotated[SetupJobFeature, Depends(FoldingDependencies.setup_job)],
) -> JobResponse:
    return await feature.handle(request=request)
