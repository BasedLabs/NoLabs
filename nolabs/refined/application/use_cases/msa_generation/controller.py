__all__ = [
    'router',
]

from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.refined.application.use_cases.msa_generation.api_models import GetJobStatusResponse, JobResponse, SetupJobRequest
from nolabs.refined.application.use_cases.msa_generation.di import MsaGenerationDependencies
from nolabs.refined.application.use_cases.msa_generation.use_cases import GetJobStatusFeature, RunJobFeature, GetJobFeature, \
    SetupJobFeature

router = APIRouter(
    prefix='/api/v1/msa-generation',
    tags=['Generate msa'],

)


@router.post('/jobs/run/{job_id}')
async def run_job(
        feature: Annotated[RunJobFeature, Depends(MsaGenerationDependencies.run_job)],
        job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}',
            summary='Get job')
async def get_job(job_id: UUID, feature: Annotated[
    GetJobFeature, Depends(MsaGenerationDependencies.get_job)]) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}/status',
            summary='Get job execution status')
async def get_job_status(job_id: UUID, feature: Annotated[
    GetJobStatusFeature, Depends(MsaGenerationDependencies.get_job_status)]) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)


@router.post('/jobs',
            summary='Setup job')
async def setup_job(request: SetupJobRequest, feature: Annotated[
    SetupJobFeature, Depends(MsaGenerationDependencies.setup_job)]) -> JobResponse:
    return await feature.handle(request=request)
