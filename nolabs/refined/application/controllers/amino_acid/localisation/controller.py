__all__ = [
    'router',
    'start',
    'jobs_metadata',
    'job',
    'delete_job',
    'update'
]


from typing import List
from uuid import UUID

from dependency_injector.wiring import Provide, inject
from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.controllers.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.controllers.amino_acid.localisation.api_models import RunJobResponse, GetJobResponse, \
    GetJobMetadataResponse, UpdateJobRequest
from nolabs.refined.application.features.amino_acid.localisation import RunLocalisationFeature, \
    GetJobFeature, DeleteJobFeature, UpdateJobFeature, GetJobsMetadataFeature

router = APIRouter(
    prefix='/api/v1/localisation',
    tags=['localisation'],

)


@router.post('/jobs/start',
             summary='Start localisation job and get probabilities of localisation of certain amino acids in cell')
@inject
async def start(
        job_id: UUID = Form(None),
        experiment_id: UUID = Form(None),
        fastas: List[UploadFile] = File(default_factory=list),
        feature: RunLocalisationFeature = Depends(Provide[RunLocalisationFeature]),
) -> RunJobResponse:
    return await feature.handle(RunAminoAcidRequest(
        job_id=job_id,
        experiment_id=experiment_id,
        fastas=fastas
    ))


@router.get('/jobs/metadata',
            summary='Get all jobs metadata by experiment')
@inject
async def jobs_metadata(experiment_id: UUID,
                        feature: GetJobsMetadataFeature = Depends(Provide[GetJobsMetadataFeature])) -> List[
    GetJobMetadataResponse]:
    return await feature.handle(experiment_id=experiment_id)


@router.get('/jobs/{job_id}',
            summary='Get job execution result')
@inject
async def job(job_id: UUID, feature: GetJobFeature = Depends(Provide[GetJobFeature])) -> GetJobResponse:
    return await feature.handle(job_id=job_id)


@router.delete('/jobs/{jod_id}',
               summary='Delete job')
@inject
async def delete_job(job_id: UUID, feature: DeleteJobFeature = Depends(Provide[DeleteJobFeature])):
    return feature.handle(job_id=job_id)


@router.patch('/jobs/{job_id}',
              summary='Update job')
@inject
async def update(job_id: UUID,
                 request: UpdateJobRequest,
                 feature: UpdateJobFeature = Depends(Provide[UpdateJobFeature])):
    return feature.handle(job_id=job_id, request=request)
