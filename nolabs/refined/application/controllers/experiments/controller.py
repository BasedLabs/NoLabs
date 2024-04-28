from typing import Annotated, List
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.refined.application.controllers.experiments.api_models import ExperimentMetadataResponse, \
    ChangeExperimentNameRequest
from nolabs.refined.application.features.experiments import GetExperimentsMetadataFeature, DeleteExperimentFeature, \
    ChangeExperimentNameFeature, CreateExperimentFeature


router = APIRouter(
    prefix='/api/v1/experiments',
    tags=['Experiments']
)


@router.get('/experiments/all')
async def experiments(
        feature: Annotated[GetExperimentsMetadataFeature, Depends(get_experiments_metadata_feature_dependency)]) -> List[ExperimentMetadataResponse]:
    return feature.handle()


@router.delete('/experiment/{experiment_id}')
async def delete_experiment(experiment_id: UUID, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]):
    return feature.handle(experiment_id)


@router.post('/experiment/name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]):
    return feature.handle(request)


@router.post('/experiment')
async def create_experiment(feature: Annotated[
    CreateExperimentFeature, Depends(create_experiment_dependency)]) -> ExperimentMetadataResponse:
    return feature.handle()
