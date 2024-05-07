__all__ = [
    'router',
]

from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.binding_pockets.api_models import JobResponse, SetupJobRequest
from nolabs.refined.application.binding_pockets.di import BindingPocketsDependencies
from nolabs.refined.application.binding_pockets.use_cases import RunJobFeature, SetupJobFeature

router = APIRouter(
    prefix='/api/v1/binding-pockets',
    tags=['Binding pockets'],

)


@router.post('/jobs/run/{job_id}',
             summary='Start binding pockets prediction job')
async def run_job(
        feature: Annotated[
            RunJobFeature, Depends(BindingPocketsDependencies.run_job)],
        job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}', summary='Get job')
async def get_job(job_id: UUID, feature: Annotated[
    RunJobFeature, Depends(BindingPocketsDependencies.get_job)]) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.post('/jobs', summary='Setup job')
async def setup_job(request: SetupJobRequest, feature: Annotated[
    SetupJobFeature, Depends(BindingPocketsDependencies.setup_job)]) -> JobResponse:
    return await feature.handle(request=request)


@router.get('/jobs/{job_id}/status', summary='Job status')
async def get_job_status(request: SetupJobRequest, feature: Annotated[
    SetupJobFeature, Depends(BindingPocketsDependencies.get_job_status)]) -> JobResponse:
    return await feature.handle(request=request)
