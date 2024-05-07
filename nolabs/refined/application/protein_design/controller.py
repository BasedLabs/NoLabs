__all__ = [
    'router',
]

from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.refined.application.protein_design.api_models import SetupJobRequest, JobResponse
from nolabs.refined.application.protein_design.di import ProteinDesignDependencies
from nolabs.refined.application.protein_design.use_cases import RunJobFeature, GetJobFeature, SetupJobFeature

router = APIRouter(
    prefix='/api/v1/protein-design',
    tags=['Protein design'],

)


@router.post('/jobs/run/{job_id}')
async def run_job(
        feature: Annotated[RunJobFeature, Depends(ProteinDesignDependencies.run_job)],
        job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}',
            summary='Get job')
async def get_job(job_id: UUID, feature: Annotated[
    GetJobFeature, Depends(ProteinDesignDependencies.get_job)]) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.post('/jobs',
            summary='Setup job')
async def get_job_status(request: SetupJobRequest, feature: Annotated[
    SetupJobFeature, Depends(ProteinDesignDependencies.setup_job)]) -> JobResponse:
    return await feature.handle(request=request)
