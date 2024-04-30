__all__ = [
    'router',
]

from typing import List, Annotated, Optional
from uuid import UUID

from dependency_injector.wiring import inject
from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.controllers.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.controllers.amino_acid.localisation.api_models import RunJobResponse, GetJobResponse, \
    GetJobMetadataResponse, UpdateJobRequest
from nolabs.refined.application.controllers.amino_acid.localisation.di import LocalisationDependencies
from nolabs.refined.application.features.amino_acid.localisation.use_cases import RunLocalisationFeature, \
    GetJobFeature, DeleteJobFeature, GetJobsMetadataFeature, UpdateJobFeature

router = APIRouter(
    prefix='/api/v1/localisation',
    tags=['localisation'],

)


@router.post('/jobs/start',
             summary='Start localisation job and get probabilities of localisation of certain amino acids in cell')
@inject
async def start(
        feature: Annotated[
            RunLocalisationFeature, Depends(LocalisationDependencies.start)],
        job_id: Optional[UUID] = Form(None),
        experiment_id: UUID = Form(),
        fastas: List[UploadFile] = File(default_factory=list)
) -> RunJobResponse:
    return await feature.handle(RunAminoAcidRequest(
        job_id=job_id,
        experiment_id=experiment_id,
        fastas=fastas
    ))


@router.get('/jobs/metadata',
            summary='Get all jobs metadata by experiment')
async def jobs_metadata(experiment_id: UUID,
                        feature: Annotated[
                            GetJobsMetadataFeature, Depends(LocalisationDependencies.job_metadata)]) -> \
        List[
            GetJobMetadataResponse]:
    return await feature.handle(experiment_id=experiment_id)


@router.get('/jobs/{job_id}',
            summary='Get job execution result')
@inject
async def job(job_id: UUID, feature: Annotated[
    GetJobFeature, Depends(LocalisationDependencies.get_job)]) -> GetJobResponse:
    return await feature.handle(job_id=job_id)


@router.delete('/jobs/{jod_id}',
               summary='Delete job')
@inject
async def delete_job(job_id: UUID,
                     feature: Annotated[DeleteJobFeature, Depends(LocalisationDependencies.delete_job)]):
    return await feature.handle(job_id=job_id)


@router.patch('/jobs/{job_id}',
             summary='Update job')
@inject
async def update(job_id: UUID,
                request: UpdateJobRequest,
                feature: Annotated[UpdateJobFeature, Depends(LocalisationDependencies.update_job)]):
   return await feature.handle(job_id=job_id, request=request)
