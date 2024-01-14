from typing import Annotated, Optional, List, Union

from fastapi import APIRouter, Depends, UploadFile, File, Form

from nolabs.api_models.solubility import RunSolubilityRequest, RunSolubilityResponse, ExperimentMetadataResponse, \
    GetExperimentResponse, ChangeExperimentNameRequest, GenerateUuidResponse
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
async def inference(
        feature: Annotated[RunSolubilityFeature, Depends(run_solubility_feature_dependency)],
        experiment_name: str = Form(),
        experiment_id: str = Form(None),
        amino_acid_sequence: str = Form(None),
        fastas: List[UploadFile] = File(default_factory=list)
) -> RunSolubilityResponse:
        return await feature.handle(RunSolubilityRequest(
            experiment_name=experiment_name,
            experiment_id=experiment_id,
            amino_acid_sequence=amino_acid_sequence,
            fastas=fastas
        ))


@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature, Depends(get_experiments_feature_dependency)]) -> List[ExperimentMetadataResponse]:
    return feature.handle()


@router.get('/get-experiment')
async def get_experiment(experiment_id: str, feature: Annotated[
    GetExperimentFeature, Depends(get_experiment_feature_dependency)]) -> GetExperimentResponse:
    return feature.handle(experiment_id)


@router.delete('/delete-experiment')
async def delete_experiment(experiment_id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]):
    return feature.handle(experiment_id)



@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]):
    return feature.handle(request)


@router.get('/generate_id')
async def generate_uuid() -> GenerateUuidResponse:
    return GenerateUuidResponse(uuid=uuid_utils.generate_uuid())
