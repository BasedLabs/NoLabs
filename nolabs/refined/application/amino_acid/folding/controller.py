__all__ = [
    'router',
]

from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.refined.application.amino_acid.folding.api_models import (JobResponse, GetJobStatusResponse,
                                                                      SetupJobRequest)
from nolabs.refined.application.amino_acid.folding.di import FoldingDependencies
from nolabs.refined.application.amino_acid.folding.use_cases import RunJobFeature, GetJobFeature, GetJobStatusFeature, \
    SetupJobFeature

router = APIRouter(
    prefix='/api/v1/folding',
    tags=['Folding'],

)


@router.post('/jobs/start/{job_id}',
             summary='Start folding job')
async def start_job(
        feature: Annotated[
            RunJobFeature, Depends(FoldingDependencies.run_job)],
        job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id)


@router.get('/jobs/{job_id}',
            summary='Get job')
async def get_job(job_id: UUID, feature: Annotated[
    GetJobFeature, Depends(FoldingDependencies.get_job)]) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}/status',
            summary='Get job execution status')
async def get_job_status(job_id: UUID, feature: Annotated[
    GetJobStatusFeature, Depends(FoldingDependencies.get_job_status)]) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)


@router.post('/jobs',
            summary='Setup job')
async def setup_job(request: SetupJobRequest, feature: Annotated[
    SetupJobFeature, Depends(FoldingDependencies.setup_job)]) -> JobResponse:
    return await feature.handle(request=request)