__all__ = [
    'router',
]

from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.use_cases.folding.use_cases import RunJobFeature
from nolabs.application.use_cases.localisation.api_models import JobResponse, SetupJobRequest, GetJobStatusResponse
from nolabs.application.use_cases.localisation.di import LocalisationDependencies
from nolabs.application.use_cases.localisation.use_cases import GetJobFeature, SetupJobFeature, GetJobStatusFeature

router = APIRouter(
    prefix='/api/v1/localisation',
    tags=['Localisation'],

)


@router.post('/jobs/run/{job_id}',
             summary='Start localisation job')
async def run_job(
        feature: Annotated[
            RunJobFeature, Depends(LocalisationDependencies.run_job)],
        job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}',
            summary='Get job')
async def get_job(job_id: UUID, feature: Annotated[
    GetJobFeature, Depends(LocalisationDependencies.get_job)]) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.post('/jobs',
            summary='Setup job')
async def setup_job(request: SetupJobRequest, feature: Annotated[
    SetupJobFeature, Depends(LocalisationDependencies.setup_job)]) -> JobResponse:
    return await feature.handle(request=request)


@router.get('/jobs/{job_id}/status',
            summary='Get job execution status')
async def get_job_status(job_id: UUID, feature: Annotated[
    GetJobStatusFeature, Depends(LocalisationDependencies.get_job_status)]) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)