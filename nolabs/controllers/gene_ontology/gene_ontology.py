from typing import Annotated, Dict

from fastapi import APIRouter, Depends

from nolabs.controllers.gene_ontology.dependencies import change_experiment_name_dependency, \
    delete_experiment_feature_dependency, get_experiment_feature_dependency, get_experiments_feature_dependency, \
    run_gene_ontology_feature_dependency
from nolabs.api_models.gene_ontology import RunGeneOntologyRequest, RunGeneOntologyResponse, ExperimentMetadataResponse, \
    GetExperimentResponse, ChangeExperimentNameRequest, GenerateUuidResponse
from nolabs.features.gene_ontology import DeleteExperimentFeature, RunGeneOntologyFeature, \
    GetExperimentsFeature, GetExperimentFeature, ChangeExperimentNameFeature
from nolabs.utils import uuid_utils

router = APIRouter(
    prefix='/api/v1/gene-ontology',
    tags=['gene-ontology']
)


@router.post('/inference')
async def inference(request: RunGeneOntologyRequest,
                    feature: Annotated[RunGeneOntologyFeature, Depends(run_gene_ontology_feature_dependency)]
                    ) -> RunGeneOntologyResponse:
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
