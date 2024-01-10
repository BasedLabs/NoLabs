from typing import Annotated, Dict

from fastapi import WebSocket, APIRouter, Depends

from nolabs.api_models.conformations import RunSimulationsRequest, RunSimulationsResponse, ExperimentMetadataResponse, \
    GetExperimentResponse, ChangeExperimentNameRequest, GenerateUuidResponse
from nolabs.controllers.common_dependencies import generate_uuid_dependency
from nolabs.controllers.conformations.dependencies import run_simulations_feature_dependency, \
    get_experiment_feature_dependency, delete_experiment_feature_dependency, \
    change_experiment_name_dependency, get_experiments_feature_dependency
from nolabs.features.conformations import DeleteExperimentFeature, RunSimulationsFeature, \
    GetExperimentsFeature, GetExperimentFeature, ChangeExperimentNameFeature
from nolabs.infrastructure.websockets import ConformationsWebsocket
from nolabs.utils import uuid_utils

logs_websocket: WebSocket | None


def conformations_websocket_dependency() -> ConformationsWebsocket:
    assert logs_websocket

    return ConformationsWebsocket(logs_websocket)


router = APIRouter(
    prefix='/api/v1/conformations',
    tags=['conformations']
)


@router.websocket("/logs")
async def ws(websocket: WebSocket):
    global logs_websocket
    await websocket.accept()
    logs_websocket = websocket


@router.post('/inference')
async def inference(request: RunSimulationsRequest,
                    feature: Annotated[RunSimulationsFeature, Depends(run_simulations_feature_dependency)]
                    ) -> RunSimulationsResponse:
    return await feature.handle(request)


@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature, Depends(get_experiments_feature_dependency)]) -> Dict[
    str, ExperimentMetadataResponse]:
    return feature.handle()


@router.get('/load-experiment')
async def get_experiment(experiment_id: str, feature: Annotated[
    GetExperimentFeature, Depends(get_experiment_feature_dependency)]) -> GetExperimentResponse:
    return feature.handle(experiment_id)


@router.delete('/delete-experiment')
async def delete_experiment(experiment_id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]) -> GetExperimentResponse:
    return feature.handle(experiment_id)


@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]) -> GetExperimentResponse:
    return feature.handle(request)


@router.get('/generate_id')
async def generate_uuid(feature: Annotated[uuid_utils.UuidUtils
, Depends(generate_uuid_dependency)]) -> GenerateUuidResponse:
    return GenerateUuidResponse(
        uuid=feature.generate_uuid()
    )
