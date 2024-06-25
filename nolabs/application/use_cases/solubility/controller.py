__all__ = [
    'router',
]

from typing import Annotated
from uuid import UUID

import fastapi

from nolabs.application.use_cases.solubility.api_models import JobResponse, SetupJobRequest, GetJobStatusResponse
from nolabs.application.use_cases.solubility.di import SolubilityDependencies
from nolabs.application.use_cases.solubility.use_cases import RunJobFeature, GetJobFeature, SetupJobFeature, \
    GetJobStatusFeature

router = fastapi.APIRouter(
    prefix='/api/v1/solubility',
    tags=['Solubility'],

)


@router.post('/jobs/run/{job_id}',
             summary='Start solubility probability determination job and get probability of protein being soluble')
async def run_job(
        feature: Annotated[
            RunJobFeature, fastapi.Depends(SolubilityDependencies.run_job)],
        job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}',
            summary='Get job')
async def get_job(job_id: UUID, feature: Annotated[
    GetJobFeature, fastapi.Depends(SolubilityDependencies.get_job)]) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.post('/jobs',
            summary='Setup job')
async def setup_job(request: SetupJobRequest, feature: Annotated[
    SetupJobFeature, fastapi.Depends(SolubilityDependencies.setup_job)]) -> JobResponse:
    return await feature.handle(request=request)


@router.get('/jobs/{job_id}/status',
            summary='Get job execution status')
async def get_job_status(job_id: UUID, feature: Annotated[
    GetJobStatusFeature, fastapi.Depends(SolubilityDependencies.get_job_status)]) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)