from typing import Annotated, Dict, List

from fastapi import APIRouter, Depends, Form, UploadFile, File

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
async def inference(feature: Annotated[RunGeneOntologyFeature, Depends(run_gene_ontology_feature_dependency)],
                    experiment_name: str = Form(),
                    experiment_id: str = Form(None),
                    amino_acid_sequence: str = Form(None),
                    fastas: List[UploadFile] = File(default_factory=list)
                    ) -> RunGeneOntologyResponse:
    return await feature.handle(RunGeneOntologyRequest(
        experiment_name=experiment_name,
        experiment_id=experiment_id,
        amino_acid_sequence=amino_acid_sequence,
        fastas=fastas
    ))


@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature, Depends(get_experiments_feature_dependency)]) -> List[
    ExperimentMetadataResponse]:
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
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]):
    return feature.handle(request)


@router.get('/generate_id')
async def generate_uuid() -> GenerateUuidResponse:
    return GenerateUuidResponse(uuid=uuid_utils.generate_uuid())
