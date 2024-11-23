__all__ = [
    'router',
]

from uuid import UUID

from fastapi import APIRouter

from nolabs.application.blast.api_models import (JobResponse, SetupJobRequest)
from nolabs.application.blast.use_cases import GetJobFeature, SetupJobFeature

router = APIRouter(
    prefix='/api/v1/blast',
    tags=['Blast']
)


@router.get('/jobs/{job_id}',
            summary='Get job')
async def get_job(job_id: UUID) -> JobResponse:
    return await GetJobFeature().handle(job_id=job_id)


@router.post('/jobs',
            summary='Setup job')
async def setup_job(request: SetupJobRequest) -> JobResponse:
    return await SetupJobFeature().handle(request=request)