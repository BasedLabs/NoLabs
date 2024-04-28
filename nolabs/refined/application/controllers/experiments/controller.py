from typing import Annotated, List
from uuid import UUID

from dependency_injector.wiring import Provide
from fastapi import APIRouter, Depends

from nolabs.refined.application.controllers.experiments.api_models import ExperimentMetadataResponse, \
    ChangeExperimentNameRequest
from nolabs.refined.application.features.experiments import GetExperimentsMetadataFeature, DeleteExperimentFeature, \
    ChangeExperimentNameFeature, CreateExperimentFeature

router = APIRouter(
    prefix='/api/v1/experiments',
    tags=['Experiments']
)


@router.get('/experiments/all',
            summary='Get all experiments'
            )
async def experiments(
        feature: Depends(Provide[GetExperimentsMetadataFeature])) -> List[ExperimentMetadataResponse]:
    return feature.handle()


@router.delete('/experiments/{experiment_id}',
               summary='Delete experiment')
async def delete_experiment(experiment_id: UUID,
                            feature: DeleteExperimentFeature = Depends(Provide[DeleteExperimentFeature])):
    return feature.handle(experiment_id)


@router.post('/experiments/name',
             summary='Change name of experiment'
             )
async def change_name(request: ChangeExperimentNameRequest,
                                 feature: ChangeExperimentNameFeature = Depends(Provide[ChangeExperimentNameFeature])):
    return feature.handle(request)


@router.post('/experiments',
             summary='Create experiment')
async def create_experiment(feature: CreateExperimentFeature = Depends(Provide[CreateExperimentFeature])) -> ExperimentMetadataResponse:
    return feature.handle()
