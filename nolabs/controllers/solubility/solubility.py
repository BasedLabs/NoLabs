from typing import Annotated, Dict

from fastapi import APIRouter, Depends

from nolabs.api_models.solubility import RunSolubilityRequest, RunSolubilityResponse, ExperimentMetadataResponse, \
    GetExperimentResponse, ChangeExperimentNameRequest, GenerateUuidResponse
from nolabs.controllers.common_dependencies import generate_uuid_dependency
from nolabs.controllers.solubility.dependencies import run_solubility_feature_dependency, \
    get_experiments_feature_dependency, \
    get_experiment_feature_dependency, delete_experiment_feature_dependency, change_experiment_name_dependency
from nolabs.features.solubility import DeleteExperimentFeature, RunSolubilityFeature, \
    GetExperimentsFeature, GetExperimentFeature, ChangeExperimentNameFeature
from nolabs.utils import uuid_utils

router = APIRouter(
    prefix='/api/v1/solubility',
    tags=['solubility']
)


@router.post('/inference')
async def inference(request: RunSolubilityRequest,
                    feature: Annotated[RunSolubilityFeature, Depends(run_solubility_feature_dependency)]
                    ) -> RunSolubilityResponse:
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
