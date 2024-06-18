__all__ = [
    'router',
]

from typing import List, Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.use_cases.jobs.api_models import GetJobMetadataResponse, UpdateJobRequest
from nolabs.application.use_cases.jobs.di import JobDependencies
from nolabs.application.use_cases.jobs.use_cases import GetJobsMetadataFeature, DeleteJobFeature, \
    UpdateJobFeature

router = APIRouter(
    prefix='/api/v1/jobs',
    tags=['Jobs - a common controller for jobs management'],

)


@router.get('/jobs/metadata',
            summary='Get all jobs metadata by experiment')
async def jobs_metadata(experiment_id: UUID,
                        feature: Annotated[
                            GetJobsMetadataFeature, Depends(JobDependencies.job_metadata)]) -> \
        List[
            GetJobMetadataResponse]:
    return await feature.handle(experiment_id=experiment_id)


@router.delete('/jobs/{jod_id}',
               summary='Delete job')
async def delete_job(job_id: UUID,
                     feature: Annotated[DeleteJobFeature, Depends(JobDependencies.delete_job)]):
    return await feature.handle(job_id=job_id)


@router.patch('/jobs/{job_id}',
              summary='Update job')
async def update(job_id: UUID,
                 request: UpdateJobRequest,
                 feature: Annotated[UpdateJobFeature, Depends(JobDependencies.update_job)]):
    return await feature.handle(job_id=job_id, request=request)
