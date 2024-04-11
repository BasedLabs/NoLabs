from typing import List

from fastapi import WebSocket, APIRouter, UploadFile, File

from nolabs.api_models.experiment import ExperimentMetadataResponse, ChangeExperimentNameRequest
from nolabs.api_models.small_molecules_design import *
from nolabs.controllers.small_molecules_design.dependencies import *
from nolabs.modules.experiment.get_experiments import GetExperimentsFeature
from nolabs.modules.small_molecules_design.get_experiment import GetExperimentFeature

logs_websocket: WebSocket | None

router = APIRouter(
    prefix='/api/v1/small-molecules-design',
    tags=['small-molecules-design']
)


@router.get('/experiment/{experiment_id}/status')
async def status(
        feature: Annotated[GetExperimentStatusFeature, Depends(experiment_status_dependency)],
        experiment_id: str
) -> GetExperimentStatusResponse:
    return await feature.handle(experiment_id=experiment_id)


@router.post("/experiment/{experiment_id}/learning")
async def learning(
        feature: Annotated[StartLearningExperimentFeature, Depends(start_learning_experiment_dependency)],
        experiment_id: str):
    return await feature.handle(experiment_id=experiment_id)


@router.post("/experiment/{experiment_id}/sampling")
async def sampling(
        feature: Annotated[StartSamplingExperimentFeature, Depends(start_sampling_experiment_dependency)],
        experiment_id: str,
        request: SamplingSizeRequest):
    return await feature.handle(experiment_id=experiment_id, request=request)


@router.get('/experiment/{experiment_id}')
async def get_experiment(feature: Annotated[GetExperimentFeature, Depends(get_experiment_dependency)],
                         experiment_id: str) -> GetExperimentResponse:
    return await feature.handle(experiment_id)


@router.post('/experiment/{experiment_id}/props')
async def save_properties(
        feature: Annotated[SavePropertiesFeature, Depends(save_properties_dependency)],
        experiment_id: str,
        center_x: float,
        center_y: float, center_z: float, size_x: float,
        size_y: float, size_z: float,
        epochs: int = 50,
        batch_size: int = 128,
        minscore: float = 0.4,
        pdb_file: UploadFile = File()):
    await feature.handle(experiment_id=experiment_id, request=ExperimentPropertiesRequest(
        pdb_file=pdb_file,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        batch_size=batch_size,
        minscore=minscore,
        epochs=epochs
    ))


@router.post("/experiment/{experiment_id}/stop")
async def stop(
        feature: Annotated[StopExperimentFeature, Depends(stop_experiment_dependency)],
        experiment_id: str):
    return await feature.handle(experiment_id)


@router.delete("/experiment/{experiment_id}")
async def delete(
        feature: Annotated[DeleteExperimentFeature, Depends(delete_experiment_dependency)],
        experiment_id: str):
    return await feature.handle(experiment_id=experiment_id)


@router.get('/experiment/{experiment_id}/logs')
async def logs(feature: Annotated[GetLogsFeature, Depends(get_logs_dependency)],
               experiment_id: str) -> LogsResponse:
    return await feature.handle(experiment_id=experiment_id)


@router.get('/experiment/{experiment_id}/smiles')
async def smiles(feature: Annotated[GetSmilesFeature, Depends(get_smiles_dependency)],
                 experiment_id: str) -> List[SmilesResponse]:
    return await feature.handle(experiment_id=experiment_id)


# Experiments

@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature, Depends(get_experiments_dependency)]) -> List[
    ExperimentMetadataResponse]:
    return feature.handle()


@router.post('/experiment/name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]):
    return feature.handle(request)


@router.post('/experiment/create')
async def create_experiment(feature: Annotated[
    CreateExperimentFeature, Depends(create_experiment_dependency)]) -> ExperimentMetadataResponse:
    return await feature.handle()
