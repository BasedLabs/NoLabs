from typing import Annotated, List

from fastapi import WebSocket, APIRouter, Depends, File, Form, UploadFile

from nolabs.api_models.experiment import ExperimentMetadataResponse, ChangeExperimentNameRequest
from nolabs.features.experiment.change_experiment_name import ChangeExperimentNameFeature
from nolabs.features.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.api_models.protein_design import RunProteinDesignRequest, RunProteinDesignResponse, \
    GetExperimentResponse
from nolabs.controllers.protein_design.dependencies import run_protein_design_feature_dependency, \
    change_experiment_name_dependency, delete_experiment_feature_dependency, get_experiment_feature_dependency, \
    get_experiments_feature_dependency, create_experiment_dependency
from nolabs.features.protein_design import RunProteinDesignFeature, GetExperimentFeature
from nolabs.features.experiment.create_experiment import CreateExperimentFeature

logs_websocket: WebSocket | None

router = APIRouter(
    prefix='/api/v1/protein-design',
    tags=['protein-design']
)


@router.post('/inference')
async def inference(feature: Annotated[RunProteinDesignFeature, Depends(run_protein_design_feature_dependency)],
                    experiment_name: str = Form(),
                    experiment_id: str = Form(None),
                    pdb_file: UploadFile = File(),
                    contig: str = Form('50'),
                    number_of_desings: int = Form(1),
                    timesteps: int = Form(None),
                    hotspots: str = Form(None),
                    ) -> RunProteinDesignResponse:
    return await feature.handle(RunProteinDesignRequest(
        experiment_name=experiment_name,
        experiment_id=experiment_id,
        pdb_file=pdb_file,
        contig=contig,
        number_of_desings=number_of_desings,
        timesteps=timesteps,
        hotspots=hotspots
    ))


@router.get('/experiments-metadata')
async def experiments(feature: Annotated[GetExperimentsMetadataFeature, Depends(get_experiments_feature_dependency)]) -> List[
    ExperimentMetadataResponse]:
    return feature.handle()


@router.get('/experiment')
async def get_experiment(id: str, feature: Annotated[
    GetExperimentFeature, Depends(get_experiment_feature_dependency)]) -> GetExperimentResponse:
    return feature.handle(id)


@router.delete('/experiment')
async def delete_experiment(id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]) -> GetExperimentResponse:
    return feature.handle(id)


@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]):
    return feature.handle(request)


@router.get('/create-experiment')
async def create_experiment(feature: Annotated[CreateExperimentFeature, Depends(create_experiment_dependency)]) -> ExperimentMetadataResponse:
    return feature.handle()
