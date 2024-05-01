__all__ = [
    'router',
]

from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.amino_acid.folding.api_models import RunJobResponse, GetJobResponse
from nolabs.refined.application.amino_acid.folding.di import FoldingDependencies
from nolabs.refined.application.amino_acid.folding.use_cases import RunFoldingFeature, \
    GetFoldingJobFeature

router = APIRouter(
    prefix='/api/v1/folding',
    tags=['Folding'],

)


@router.post('/jobs/start',
             summary='Start folding job and get pdb file result')
async def start(
        feature: Annotated[
            RunFoldingFeature, Depends(FoldingDependencies.start)],
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
    GetFoldingJobFeature, Depends(FoldingDependencies.get_job)]) -> GetJobResponse:
    return await feature.handle(job_id=job_id)
