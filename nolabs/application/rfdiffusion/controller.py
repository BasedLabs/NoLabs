__all__ = [
    'router',
]

from uuid import UUID

from fastapi import APIRouter

from nolabs.application.rfdiffusion.api_models import SetupJobRequest, JobResponse
from nolabs.application.rfdiffusion.use_cases import RunJobFeature, GetJobFeature, SetupJobFeature

router = APIRouter(
    prefix='/api/v1/rfdiffusion',
    tags=['Rfdiffusion'],

)


@router.post('/jobs/run/{job_id}')
async def run_job(job_id: UUID) -> JobResponse:
    return await RunJobFeature().handle(job_id=job_id)


@router.get('/jobs/{job_id}', summary='Get job')
async def get_job(job_id: UUID) -> JobResponse:
    return await GetJobFeature().handle(job_id=job_id)


@router.post('/jobs', summary='Setup job')
async def setup_job(request: SetupJobRequest) -> JobResponse:
    return await SetupJobFeature().handle(request=request)
