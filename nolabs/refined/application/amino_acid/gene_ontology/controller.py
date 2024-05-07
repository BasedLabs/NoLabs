__all__ = [
    'router',
]

from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.refined.application.amino_acid.gene_ontology.api_models import JobResponse, SetupJobRequest
from nolabs.refined.application.amino_acid.gene_ontology.di import GeneOntologyDependencies
from nolabs.refined.application.amino_acid.gene_ontology.use_cases import RunJobFeature, GetJobFeature, SetupJobFeature

router = APIRouter(
    prefix='/api/v1/gene-ontology',
    tags=['Gene ontology'],

)


@router.post('/jobs/run/{job_id}',
             summary='Start job')
async def start_job(
        feature: Annotated[
            RunJobFeature, Depends(GeneOntologyDependencies.run_job)],
        job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.get('/jobs/{job_id}',
            summary='Get job')
async def get_job(job_id: UUID, feature: Annotated[
    GetJobFeature, Depends(GeneOntologyDependencies.get_job)]) -> JobResponse:
    return await feature.handle(job_id=job_id)


@router.post('/jobs',
             summary='Get job')
async def setup_job(request: SetupJobRequest, feature: Annotated[
    SetupJobFeature, Depends(GeneOntologyDependencies.setup_job)]) -> JobResponse:
    return await feature.handle(request=request)
