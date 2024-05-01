__all__ = [
    'router',
]

from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.amino_acid.gene_ontology.api_models import RunJobResponse, GetJobResponse
from nolabs.refined.application.amino_acid.gene_ontology.di import GeneOntologyDependencies
from nolabs.refined.application.amino_acid.gene_ontology.use_cases import RunGeneOntologyFeature, \
    GetGeneOntologyJobFeature

router = APIRouter(
    prefix='/api/v1/gene-ontology',
    tags=['Gene ontology'],

)


@router.post('/jobs/start',
             summary='Start job and get hierarchical view of protein gene ontology')
async def start(
        feature: Annotated[
            RunGeneOntologyFeature, Depends(GeneOntologyDependencies.start)],
        job_id: Optional[UUID] = Form(None),
        experiment_id: UUID = Form(),
        fastas: List[UploadFile] = File(default_factory=list)
) -> RunJobResponse:
    return await feature.handle(RunAminoAcidRequest(
        job_id=job_id,
        experiment_id=experiment_id,
        fastas=fastas
    ))


@router.get('/jobs/{job_id}',
            summary='Get job execution result')
async def job(job_id: UUID, feature: Annotated[
    GetGeneOntologyJobFeature, Depends(GeneOntologyDependencies.get_job)]) -> GetJobResponse:
    return await feature.handle(job_id=job_id)
