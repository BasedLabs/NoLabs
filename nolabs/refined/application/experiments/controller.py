__all__ = [
    'router'
]

from typing import List, Annotated
from uuid import UUID

from fastapi import APIRouter, Depends
from nolabs.refined.application.experiments.api_models import ExperimentMetadataResponse, \
    UpdateExperimentRequest
from nolabs.refined.application.experiments.di import ExperimentsDependencies
from nolabs.refined.application.experiments.use_cases import GetExperimentsMetadataFeature, \
    DeleteExperimentFeature, UpdateExperimentFeature, CreateExperimentFeature

router = APIRouter(
    prefix='/api/v1/experiments',
    tags=['Experiments']
)


@router.get('/experiments/all',
            summary='Get all experiments'
            )
async def experiments(
        feature: Annotated[
            GetExperimentsMetadataFeature, Depends(ExperimentsDependencies.experiments)]) -> List[ExperimentMetadataResponse]:
    return feature.handle()


@router.delete('/experiments/{experiment_id}',
               summary='Delete experiment')
async def delete_experiment(experiment_id: UUID,
                            feature: Annotated[
                                DeleteExperimentFeature, Depends(ExperimentsDependencies.delete_experiment)]):
    return feature.handle(experiment_id)


@router.patch('/experiments',
              summary='Update experiment'
              )
async def update_experiment(request: UpdateExperimentRequest,
                            feature: Annotated[
                                UpdateExperimentFeature, Depends(ExperimentsDependencies.update_experiment)]):
    return feature.handle(request)


@router.post('/experiments',
             summary='Create experiment')
async def create_experiment(feature: Annotated[
    CreateExperimentFeature, Depends(ExperimentsDependencies.create_experiment)]) -> ExperimentMetadataResponse:
    return feature.handle()
