from typing import Annotated, Dict

from fastapi import WebSocket, APIRouter, Depends

from nolabs.features.protein_design import GetResultsFeature, DeleteExperimentFeature, RunProteinDesignFeature, \
    GetExperimentsFeature, GetExperimentFeature, ChangeExperimentNameFeature
from nolabs.infrastructure.settings import Settings
from nolabs.api_models.protein_design import RunProteinDesignRequest, RunProteinDesignResponse, ExperimentMetadataResponse, \
    GetExperimentResponse, ChangeExperimentNameRequest, GenerateUuidResponse
from nolabs.infrastructure.websockets import ConformationsWebsocket
from nolabs.utils import datetime_utils, uuid_utils

logs_websocket: WebSocket | None


def conformations_websocket_dependency() -> ConformationsWebsocket:
    assert logs_websocket

    return ConformationsWebsocket(logs_websocket)


router = APIRouter(
    prefix='conformations',
    tags=['conformations'],
    dependencies=[Depends(Settings),
                  Depends(datetime_utils.DateTimeUtils),
                  Depends(uuid_utils.UuidUtils),
                  Depends(RunProteinDesignFeature),
                  Depends(GetResultsFeature),
                  Depends(DeleteExperimentFeature),
                  Depends(conformations_websocket_dependency),
                  Depends(GetExperimentsFeature),
                  Depends(GetExperimentFeature),
                  Depends(ChangeExperimentNameFeature)]
)


@router.websocket("/logs")
async def ws(websocket: WebSocket):
    global logs_websocket
    await websocket.accept()
    logs_websocket = websocket


@router.post('/inference')
async def inference(request: RunProteinDesignFeature,
                    feature: Annotated[RunProteinDesignFeature, Depends(RunProteinDesignFeature)]
                    ) -> RunProteinDesignFeature:
    return await feature.handle(request)


@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature, Depends(GetExperimentsFeature)]) -> Dict[
    str, ExperimentMetadataResponse]:
    return feature.handle()


@router.get('/load-experiment')
async def get_experiment(experiment_id: str, feature: Annotated[
    GetExperimentFeature, Depends(GetExperimentFeature)]) -> GetExperimentResponse:
    return feature.handle(experiment_id)


@router.delete('/delete-experiment')
async def delete_experiment(experiment_id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(DeleteExperimentFeature)]) -> GetExperimentResponse:
    return feature.handle(experiment_id)


@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(ChangeExperimentNameFeature)]) -> GetExperimentResponse:
    return feature.handle(request)


@router.get('/generate_id')
async def generate_uuid(feature: Annotated[uuid_utils.UuidUtils
, Depends(uuid_utils.UuidUtils)]) -> GenerateUuidResponse:
    return GenerateUuidResponse(
        uuid=feature.generate_uuid()
    )
