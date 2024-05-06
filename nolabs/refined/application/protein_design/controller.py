__all__ = [
    'router',
]

from typing import List, Annotated, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Form, File, UploadFile

from nolabs.refined.application.protein_design.api_models import RunProteinDesignRequest, GetJobResponse, \
    RunProteinDesignResponse
from nolabs.refined.application.protein_design.di import ProteinDesignDependencies
from nolabs.refined.application.protein_design.use_cases import RunProteinDesignFeature, \
    GetProteinDesignJobFeature

router = APIRouter(
    prefix='/api/v1/protein-design',
    tags=['Protein design'],

)


@router.post('/jobs/start')
async def inference(feature: Annotated[RunProteinDesignFeature, Depends(ProteinDesignDependencies.start)],
                    job_id: UUID = Form(None),
                    experiment_id: UUID = Form(),
                    pdb_file: UploadFile = File(),
                    contig: str = Form('50'),
                    number_of_designs: int = Form(1),
                    timesteps: int = Form(None),
                    hotspots: str = Form(None),
                    ) -> RunProteinDesignResponse:
    return await feature.handle(RunProteinDesignRequest(
        job_id=job_id,
        experiment_id=experiment_id,
        pdb_file=pdb_file,
        contig=contig,
        number_of_designs=number_of_designs,
        timesteps=timesteps,
        hotspots=hotspots
    ))


@router.get('/jobs/{job_id}')
async def get_experiment(experiment_id: UUID, job_id: UUID, feature: Annotated[
    GetProteinDesignJobFeature, Depends(ProteinDesignDependencies.get_job)]) -> GetJobResponse:
    return await feature.handle(experiment_id=experiment_id, job_id=job_id)
