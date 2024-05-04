from typing import Annotated, Union, List
from uuid import UUID

from fastapi import APIRouter, Depends, UploadFile, File, Form

from nolabs.refined.application.conformations.api_models import IntegratorsRequest, RunSimulationsResponse, \
    RunSimulationsRequest, GetJobResponse
from nolabs.refined.application.conformations.use_cases import RunConformationsJobFeature, GetConformationsJobFeature
from nolabs.refined.application.conformations.di import ConformationsDependencies

router = APIRouter(
    prefix='/api/v1/conformations',
    tags=['conformations']
)


@router.post('/jobs/start')
async def inference(
        feature: Annotated[RunConformationsJobFeature, Depends(ConformationsDependencies.start)],
        pdb_file: UploadFile = File(),
        job_id: UUID = Form(None),
        experiment_id: UUID = Form(),
        total_frames: int = Form(10000),
        temperature_k: float = Form(273.15),
        take_frame_every: int = Form(1000),
        step_size: float = Form(0.002),
        replace_non_standard_residues: bool = Form(default=False),
        add_missing_atoms: bool = Form(default=False),
        add_missing_hydrogens: bool = Form(True),
        friction_coeff: float = Form(1.0),
        ignore_missing_atoms: bool = Form(default=False),
        integrator: IntegratorsRequest = Form(default=IntegratorsRequest.langevin),
) -> RunSimulationsResponse:
    return await feature.handle(RunSimulationsRequest(
        pdb_file=pdb_file,
        job_id=job_id,
        experiment_id=experiment_id,
        total_frames=total_frames,
        temperature_k=temperature_k,
        take_frame_every=take_frame_every,
        step_size=step_size,
        replace_non_standard_residues=replace_non_standard_residues,
        add_missing_atoms=add_missing_atoms,
        add_missing_hydrogens=add_missing_hydrogens,
        friction_coeff=friction_coeff,
        ignore_missing_atoms=ignore_missing_atoms,
        integrator=integrator
    ))


@router.get('/jobs/{job_id}')
async def get_experiment(job_id: UUID, feature: Annotated[
    GetConformationsJobFeature, Depends(ConformationsDependencies.get_job)]) -> GetJobResponse:
    return await feature.handle(job_id)