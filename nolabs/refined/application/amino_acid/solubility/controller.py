__all__ = [
    'router',
]

from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.amino_acid.solubility.api_models import RunJobResponse, GetJobResponse
from nolabs.refined.application.amino_acid.solubility.di import SolubilityDependencies
from nolabs.refined.application.amino_acid.solubility.use_cases import RunSolubilityJobFeature, \
    GetSolubilityJobFeature

router = APIRouter(
    prefix='/api/v1/localisation',
    tags=['Localisation'],

)


@router.post('/jobs/start',
             summary='Start solubility probability determination job and get probability of protein being soluble')
async def start(
        feature: Annotated[
            RunSolubilityJobFeature, Depends(SolubilityDependencies.start)],
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
    GetSolubilityJobFeature, Depends(SolubilityDependencies.get_job)]) -> GetJobResponse:
    return await feature.handle(job_id=job_id)
