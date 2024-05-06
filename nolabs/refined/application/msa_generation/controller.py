__all__ = [
    'router',
]

from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.msa_generation.api_models import PredictMsaRequest, PredictMsaResponse, \
    GetJobStatusResponse
from nolabs.refined.application.msa_generation.di import MsaGenerationDependencies
from nolabs.refined.application.msa_generation.use_cases import RunMsaGenerationJobFeature, GetJobStatusFeature, \
    GetMsaPredictionJobResultFeature

router = APIRouter(
    prefix='/api/v1/msa-generation',
    tags=['Generate msa'],

)


@router.post('/jobs/start',
             summary='Start msa generation job and get msa result')
async def start(
        feature: Annotated[
            RunMsaGenerationJobFeature, Depends(MsaGenerationDependencies.start)],
        request: PredictMsaRequest
) -> PredictMsaResponse:
    return await feature.handle(request=request)


@router.get('/jobs/{job_id}',
            summary='Get job execution result')
async def job(job_id: UUID, feature: Annotated[
    GetMsaPredictionJobResultFeature, Depends(MsaGenerationDependencies.get_job)]) -> PredictMsaResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}/status',
            summary='Get job execution status')
async def job(job_id: UUID, feature: Annotated[
    GetJobStatusFeature, Depends(MsaGenerationDependencies.get_status)]) -> GetJobStatusResponse:
    return await feature.handle(job_id=job_id)
