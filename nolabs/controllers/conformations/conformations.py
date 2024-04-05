from typing import Annotated, Union, List

from fastapi import APIRouter, Depends, UploadFile, File, Form

from nolabs.api_models.conformations import RunSimulationsRequest, RunSimulationsResponse, \
    GetExperimentResponse, IntegratorsRequest
from nolabs.api_models.experiment import ChangeExperimentNameRequest, ExperimentMetadataResponse
from nolabs.api_models.problem_details import ProblemDetailsResponse
from nolabs.controllers.conformations.dependencies import run_simulations_feature_dependency, \
    get_experiment_feature_dependency, delete_experiment_feature_dependency, \
    change_experiment_name_dependency, get_experiments_feature_dependency, create_experiment_dependency
from nolabs.modules.conformations.get_experiment import GetExperimentFeature
from nolabs.modules.conformations.run_simulations import RunSimulationsFeature
from nolabs.modules.experiment.change_experiment_name import ChangeExperimentNameFeature
from nolabs.modules.experiment.create_experiment import CreateExperimentFeature
from nolabs.modules.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.modules.experiment.get_experiments import GetExperimentsFeature

router = APIRouter(
    prefix='/api/v1/conformations',
    tags=['conformations']
)


@router.post('/inference')
async def inference(
        feature: Annotated[RunSimulationsFeature, Depends(run_simulations_feature_dependency)],
        pdb_file: UploadFile = File(),
        experiment_name: str = Form(),
        experiment_id: str = Form(None),
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
    global logs_websocket

    return await feature.handle(RunSimulationsRequest(
        pdb_file=pdb_file,
        experiment_name=experiment_name,
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
        integrator=integrator,
    ))


@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature, Depends(get_experiments_feature_dependency)]) -> List[
    ExperimentMetadataResponse]:
    return feature.handle()


@router.get('/get-experiment')
async def get_experiment(experiment_id: str, feature: Annotated[
    GetExperimentFeature, Depends(get_experiment_feature_dependency)]) -> GetExperimentResponse:
    return await feature.handle(experiment_id)


@router.delete('/delete-experiment')
async def delete_experiment(experiment_id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]):
    return feature.handle(experiment_id)


@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]):
    return feature.handle(request)

@router.get('/create-experiment')
async def create_experiment(feature: Annotated[CreateExperimentFeature, Depends(create_experiment_dependency)]) -> ExperimentMetadataResponse:
    return await feature.handle()