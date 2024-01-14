from typing import Annotated, Dict

from fastapi import WebSocket, APIRouter, Depends

from nolabs.api_models.protein_design import RunProteinDesignRequest, RunProteinDesignResponse, \
    ExperimentMetadataResponse, \
    GetExperimentResponse, ChangeExperimentNameRequest, GenerateUuidResponse
from nolabs.controllers.protein_design.dependencies import run_protein_design_feature_dependency, \
    change_experiment_name_dependency, delete_experiment_feature_dependency, get_experiment_feature_dependency, \
    get_experiments_feature_dependency
from nolabs.features.protein_design import DeleteExperimentFeature, RunProteinDesignFeature, \
    GetExperimentsFeature, GetExperimentFeature, ChangeExperimentNameFeature
from nolabs.utils import uuid_utils

logs_websocket: WebSocket | None


router = APIRouter(
    prefix='/api/v1/conformations',
    tags=['conformations']
)


@router.post('/inference')
async def inference(request: RunProteinDesignRequest,
                    feature: Annotated[RunProteinDesignFeature, Depends(run_protein_design_feature_dependency)]
                    ) -> RunProteinDesignResponse:
    return await feature.handle(request)


@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature, Depends(get_experiments_feature_dependency)]) -> Dict[
    str, ExperimentMetadataResponse]:
    return feature.handle()


@router.get('/get-experiment')
async def get_experiment(experiment_id: str, feature: Annotated[
    GetExperimentFeature, Depends(get_experiment_feature_dependency)]) -> GetExperimentResponse:
    return feature.handle(experiment_id)


@router.delete('/delete-experiment')
@router.delete('/delete-experiment')
async def delete_experiment(experiment_id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]) -> GetExperimentResponse:
    return feature.handle(experiment_id)


@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]) -> GetExperimentResponse:
    return feature.handle(request)


@router.get('/generate_id')
async def generate_uuid() -> GenerateUuidResponse:
    return GenerateUuidResponse(uuid=uuid_utils.generate_uuid())
