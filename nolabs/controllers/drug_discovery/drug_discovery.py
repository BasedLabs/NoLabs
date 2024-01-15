from typing import List

from fastapi import APIRouter, Depends, UploadFile, File, Form
from typing import Annotated
from nolabs.controllers.drug_discovery.dependencies import (
    get_experiments_feature_dependency,
    add_experiment_feature_dependency,
    delete_experiment_feature_dependency,
    change_experiment_name_dependency,
    upload_target_dependency,
    delete_target_dependency,
    get_targets_list_dependency,
    upload_ligand_dependency,
    delete_ligand_dependency,
    get_ligands_list_dependency,
    get_binding_pocket_dependency,
    predict_binding_pocket_dependency,
    predict_folding_dependency,
    get_folded_structure_dependency,
    generate_msa_dependency,
    predict_docking_dependency
)
from nolabs.features.drug_discovery import (
    GetExperimentsFeature,
    DeleteExperimentFeature,
    ChangeExperimentNameFeature,
    AddExperimentFeature
)
from nolabs.features.drug_discovery.target_management import UploadTargetFeature, DeleteTargetFeature, \
    GetTargetsListFeature
from nolabs.features.drug_discovery.get_binding_pocket import GetBindingPocketFeature
from nolabs.features.drug_discovery.predict_binding_pocket import PredictBindingPocketFeature
from nolabs.features.drug_discovery.predict_light_folding import PredictFoldingFeature
from nolabs.features.drug_discovery.get_folding import GetFoldedStructureFeature
from nolabs.features.drug_discovery.generate_msa import GenerateMsaFeature
from nolabs.features.drug_discovery.ligand_management import UploadLigandFeature, DeleteLigandFeature, \
    GetLigandsListFeature
from nolabs.features.drug_discovery.predict_docking import PredictDockingFeature

from nolabs.api_models.drug_discovery import (
    ChangeExperimentNameRequest,
    UploadTargetRequest,
    DeleteTargetRequest,
    GetTargetsListRequest,
    UploadLigandRequest,
    DeleteLigandRequest,
    GetLigandsListRequest,
    GetTargetBindingPocketRequest,
    PredictBindingPocketRequest,
    PredictMsaRequest,
    GetFoldingRequest,
    PredictFoldingRequest,
    DockingRequest, ExperimentMetadataResponse
)

router = APIRouter(
    prefix='/api/v1/drug_discovery',
    tags=['drug_discovery']
)


@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature,
                      Depends(get_experiments_feature_dependency)]) -> List[ExperimentMetadataResponse]:
    return feature.handle()

@router.post('/add-experiment')
async def add_experiment(feature: Annotated[
    AddExperimentFeature, Depends(add_experiment_feature_dependency)]) -> ExperimentMetadataResponse:
    return feature.handle()

@router.delete('/delete-experiment')
async def delete_experiment(experiment_id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]):
    return feature.handle(experiment_id)


@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]):
    return feature.handle(request)


# Target Management
@router.post('/upload-target')
async def upload_target(feature: Annotated[UploadTargetFeature, Depends(upload_target_dependency)],
                        experiment_id: str = Form(),
                        fasta: UploadFile = File()
                        ):
    return feature.handle(UploadTargetRequest(experiment_id, fasta))


@router.delete('/delete-target')
async def delete_target(request: DeleteTargetRequest, feature: Annotated[
    DeleteTargetFeature, Depends(delete_target_dependency)]):
    return feature.handle(request)


@router.get('/get-targets-list')
async def get_targets_list(experiment_id: str, feature: Annotated[
    GetTargetsListFeature, Depends(get_targets_list_dependency)]):
    return feature.handle(GetTargetsListRequest(experiment_id))


# Ligand Management
@router.post('/upload-ligand')
async def upload_ligand(feature: Annotated[UploadLigandFeature, Depends(upload_ligand_dependency)],
                        experiment_id: str = Form(),
                        target_id: str = Form(),
                        sdf_file: UploadFile = File()
                        ):
    return feature.handle(UploadLigandRequest(experiment_id, target_id, sdf_file))


@router.delete('/delete-ligand')
async def delete_ligand(feature: Annotated[
                        DeleteLigandFeature, Depends(delete_ligand_dependency)],
                        experiment_id: str,
                        target_id: str,
                        ligand_id: str):
    return feature.handle(DeleteLigandRequest(experiment_id, target_id, ligand_id))


@router.get('/get-ligands-list')
async def get_ligands_list(feature: Annotated[
    GetLigandsListFeature, Depends(get_ligands_list_dependency)],
                           experiment_id: str,
                           target_id: str):
    return feature.handle(GetLigandsListRequest(experiment_id, target_id))


# Binding Pocket Management
@router.get('/get-target-binding-pocket')
async def get_target_binding_pocket(feature: Annotated[GetBindingPocketFeature, Depends(get_binding_pocket_dependency)],
                                    experiment_id: str,
                                    target_id: str
                                    ):
    return feature.handle(GetTargetBindingPocketRequest(experiment_id, target_id))


@router.post('/predict-binding-pocket')
async def predict_binding_pocket(request: PredictBindingPocketRequest, feature: Annotated[
    PredictBindingPocketFeature, Depends(predict_binding_pocket_dependency)]):
    return feature.handle(request)


# Folding and MSA Prediction
@router.post('/predict-msa')
async def predict_msa(request: PredictMsaRequest, feature: Annotated[
    GenerateMsaFeature, Depends(generate_msa_dependency)]):
    return feature.handle(request)


@router.get('/get-folded-structure')
async def check_folding_exist(request: GetFoldingRequest,
                              feature: Annotated[GetFoldedStructureFeature, Depends(get_folded_structure_dependency)]):
    return feature.handle(request)


@router.post('/predict-folding')
async def predict_folding(request: PredictFoldingRequest, feature: Annotated[
    PredictFoldingFeature, Depends(predict_folding_dependency)]):
    return feature.handle(request)


# Docking
@router.post('/predict-docking')
async def perform_docking(request: DockingRequest, feature: Annotated[
    PredictDockingFeature, Depends(predict_docking_dependency)]):
    return feature.handle(request)
