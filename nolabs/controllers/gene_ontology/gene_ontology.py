from typing import Annotated, Dict, List

from fastapi import APIRouter, Depends, Form, UploadFile, File

from nolabs.controllers.gene_ontology.dependencies import change_experiment_name_dependency, \
    delete_experiment_feature_dependency, get_experiment_feature_dependency, get_experiments_feature_dependency, \
    run_gene_ontology_feature_dependency, create_experiment_dependency
from nolabs.api_models.gene_ontology import RunGeneOntologyRequest, RunGeneOntologyResponse, \
    GetExperimentResponse
from nolabs.features.experiment.create_experiment import CreateExperimentFeature
from nolabs.features.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.features.experiment.change_experiment_name import ChangeExperimentNameFeature
from nolabs.api_models.experiment import ChangeExperimentNameRequest, ExperimentMetadataResponse
from nolabs.features.experiment.get_experiments import GetExperimentsFeature
from nolabs.features.gene_ontology.get_experiment import GetExperimentFeature
from nolabs.features.gene_ontology.run_gene_ontology import RunGeneOntologyFeature
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
    return await feature.handle(experiment_id)


@router.delete('/delete-experiment')
async def delete_experiment(experiment_id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]) -> GetExperimentResponse:
    return feature.handle(experiment_id)


@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]):
    return feature.handle(request)


@router.get('/create-experiment')
async def create_experiment(feature: Annotated[CreateExperimentFeature, Depends(create_experiment_dependency)]) -> ExperimentMetadataResponse:
    return feature.handle()
