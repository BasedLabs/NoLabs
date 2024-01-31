from typing import List

from fastapi import APIRouter, Depends, UploadFile, File, Form
from typing import Annotated

from nolabs.api_models.experiment import ExperimentMetadataResponse, ChangeExperimentNameRequest
from nolabs.controllers.drug_discovery.dependencies import (
    get_experiments_feature_dependency,
    delete_experiment_feature_dependency,
    change_experiment_name_dependency,
    upload_target_dependency,
    delete_target_dependency,
    get_targets_list_dependency,
    get_target_data_dependency,
    get_ligand_data_dependency,
    upload_ligand_dependency,
    delete_ligand_dependency,
    get_ligands_list_dependency,
    get_binding_pocket_dependency,
    predict_binding_pocket_dependency,
    predict_folding_dependency,
    get_folded_structure_dependency,
    generate_msa_dependency,
    register_docking_job_dependency,
    check_msa_running_dependency,
    check_p2rank_running_dependency,
    check_umol_running_dependency,
    check_result_data_available_dependency,
    get_results_list_for_target_ligand_feature,
    get_all_results_list_dependency,
    predict_docking_dependency,
    get_docking_result_dependency, create_experiment_dependency
)
from nolabs.features.drug_discovery.result_management import CheckResultDataAvailableFeature, \
    GetAllResultsListFeature, GetResultsListForTargetLigandFeature
from nolabs.features.drug_discovery.target_management import UploadTargetFeature, DeleteTargetFeature, \
    GetTargetsListFeature, GetTargetDataFeature
from nolabs.features.drug_discovery.get_binding_pocket import GetBindingPocketFeature
from nolabs.features.drug_discovery.predict_binding_pocket import PredictBindingPocketFeature
from nolabs.features.drug_discovery.predict_light_folding import PredictFoldingFeature
from nolabs.features.drug_discovery.get_folding import GetFoldedStructureFeature
from nolabs.features.drug_discovery.generate_msa import GenerateMsaFeature
from nolabs.features.drug_discovery.ligand_management import UploadLigandFeature, DeleteLigandFeature, \
    GetLigandsListFeature, GetLigandDataFeature
from nolabs.features.drug_discovery.register_docking_job import RegisterDockingJobFeature
from nolabs.features.drug_discovery.predict_docking import PredictDockingFeature
from nolabs.features.experiment.create_experiment import CreateExperimentFeature
from nolabs.features.drug_discovery.get_results import GetDockingResultsFeature
from nolabs.features.drug_discovery.progress_management import CheckMsaRunningFeature, \
    CheckP2RankRunningFeature, CheckUmolRunningFeature
from nolabs.features.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.features.experiment.change_experiment_name import ChangeExperimentNameFeature

from nolabs.api_models.drug_discovery import (
    UploadTargetRequest,
    UploadTargetResponse,
    DeleteTargetRequest,
    DeleteTargetResponse,
    GetTargetsListRequest,
    GetTargetDataRequest,
    GetTargetDataResponse,
    UploadLigandRequest,
    UploadLigandResponse,
    DeleteLigandRequest,
    DeleteLigandResponse,
    GetLigandsListRequest,
    GetLigandDataResponse,
    GetTargetBindingPocketRequest,
    GetTargetBindingPocketResponse,
    PredictBindingPocketRequest,
    PredictBindingPocketResponse,
    PredictMsaRequest,
    PredictMsaResponse,
    GetFoldingRequest,
    GetFoldingResponse,
    PredictFoldingRequest,
    PredictFoldingResponse,
    CheckJobIsRunningRequest,
    CheckJobIsRunningResponse,
    RegisterDockingJobRequest,
    RegisterDockingJobResponse,
    RunDockingJobRequest,
    RunDockingJobResponse,
    CheckResultDataAvailableRequest,
    CheckResultDataAvailableResponse,
    GetDockingResultDataRequest,
    GetDockingResultDataResponse,
    TargetMetaData,
    LigandMetaData,
    GetLigandDataRequest,
    GetResultsListForTargetLigandRequest,
    GetResultsListForTargetLigandResponse,
    GetAllResultsListRequest,
    GetAllResultsListResponse
)
from nolabs.features.experiment.get_experiments import GetExperimentsFeature

router = APIRouter(
    prefix='/api/v1/drug_discovery',
    tags=['drug_discovery']
)


@router.get('/experiments')
async def experiments(feature: Annotated[GetExperimentsFeature,
Depends(get_experiments_feature_dependency)]) -> List[ExperimentMetadataResponse]:
    return feature.handle()


@router.get('/create-experiment')
async def create_experiment(feature: Annotated[
    CreateExperimentFeature, Depends(create_experiment_dependency)]) -> ExperimentMetadataResponse:
    return await feature.handle()


@router.delete('/delete-experiment')
async def delete_experiment(experiment_id: str, feature: Annotated[
    DeleteExperimentFeature, Depends(delete_experiment_feature_dependency)]):
    return feature.handle(experiment_id)


@router.post('/change-experiment-name')
async def change_experiment_name(request: ChangeExperimentNameRequest, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]) -> ExperimentMetadataResponse:
    return feature.handle(request)


# Target Management
@router.post('/upload-target')
async def upload_target(feature: Annotated[UploadTargetFeature, Depends(upload_target_dependency)],
                        experiment_id: str = Form(),
                        fasta: UploadFile = File()
                        ) -> UploadTargetResponse:
    return feature.handle(UploadTargetRequest(experiment_id, fasta))


@router.delete('/delete-target')
async def delete_target(feature: Annotated[
    DeleteTargetFeature, Depends(delete_target_dependency)],
                        experiment_id: str,
                        target_id: str) -> DeleteTargetResponse:
    return feature.handle(DeleteTargetRequest(experiment_id, target_id))


@router.get('/get-targets-list')
async def get_targets_list(experiment_id: str, feature: Annotated[
    GetTargetsListFeature, Depends(get_targets_list_dependency)]) -> List[TargetMetaData]:
    return feature.handle(GetTargetsListRequest(experiment_id)).targets


@router.get('/get-target-data')
async def get_targets_list(feature: Annotated[
    GetTargetDataFeature, Depends(get_target_data_dependency)],
                           experiment_id: str,
                           target_id: str,
                           ) -> GetTargetDataResponse:
    result = feature.handle(GetTargetDataRequest(experiment_id, target_id))
    return GetTargetDataResponse(protein_sequence=result.protein_sequence, protein_pdb=result.protein_pdb)


# Ligand Management
@router.post('/upload-ligand')
async def upload_ligand(feature: Annotated[UploadLigandFeature, Depends(upload_ligand_dependency)],
                        experiment_id: str = Form(),
                        target_id: str = Form(),
                        sdf_file: UploadFile = File()
                        ) -> UploadLigandResponse:
    return feature.handle(UploadLigandRequest(experiment_id, target_id, sdf_file))


@router.delete('/delete-ligand')
async def delete_ligand(feature: Annotated[
    DeleteLigandFeature, Depends(delete_ligand_dependency)],
                        experiment_id: str,
                        target_id: str,
                        ligand_id: str) -> DeleteLigandResponse:
    return feature.handle(DeleteLigandRequest(experiment_id, target_id, ligand_id))


@router.get('/get-ligand-data')
async def get_ligand_data(feature: Annotated[
    GetLigandDataFeature, Depends(get_ligand_data_dependency)],
                          experiment_id: str,
                          target_id: str,
                          ligand_id: str) -> GetLigandDataResponse:
    return feature.handle(GetLigandDataRequest(experiment_id, target_id, ligand_id))


@router.get('/get-ligands-list')
async def get_ligands_list(feature: Annotated[
    GetLigandsListFeature, Depends(get_ligands_list_dependency)],
                           experiment_id: str,
                           target_id: str) -> List[LigandMetaData]:
    return feature.handle(GetLigandsListRequest(experiment_id, target_id)).ligands


# Binding Pocket Management
@router.get('/get-target-binding-pocket')
async def get_target_binding_pocket(feature: Annotated[GetBindingPocketFeature, Depends(get_binding_pocket_dependency)],
                                    experiment_id: str,
                                    target_id: str
                                    ) -> GetTargetBindingPocketResponse:
    return feature.handle(GetTargetBindingPocketRequest(experiment_id, target_id))


@router.post('/predict-binding-pocket')
async def predict_binding_pocket(feature: Annotated[
    PredictBindingPocketFeature, Depends(predict_binding_pocket_dependency)],
                                 experiment_id: str,
                                 target_id: str
                                 ) -> PredictBindingPocketResponse:
    return feature.handle(PredictBindingPocketRequest(experiment_id, target_id))


# Folding and MSA Prediction
@router.post('/predict-msa')
async def predict_msa(feature: Annotated[
    GenerateMsaFeature, Depends(generate_msa_dependency)],
                      experiment_id: str,
                      target_id: str
                      ) -> PredictMsaResponse:
    return feature.handle(PredictMsaRequest(experiment_id, target_id))


@router.get('/get-folded-structure')
async def get_folded_structure(experiment_id: str,
                               target_id: str,
                               feature: Annotated[GetFoldedStructureFeature,
                               Depends(get_folded_structure_dependency)]) -> GetFoldingResponse:
    return feature.handle(GetFoldingRequest(experiment_id, target_id))


@router.post('/predict-folding')
async def predict_folding(feature: Annotated[
    PredictFoldingFeature, Depends(predict_folding_dependency)],
                          experiment_id: str,
                          target_id: str
                          ) -> PredictFoldingResponse:
    return feature.handle(PredictFoldingRequest(experiment_id, target_id))


@router.post('/register-docking-job')
async def register_docking_job(request: RegisterDockingJobRequest, feature: Annotated[RegisterDockingJobFeature
, Depends(register_docking_job_dependency)]) -> RegisterDockingJobResponse:
    return feature.handle(request)


@router.get('/check-msa-job-running')
async def check_msa_job_running(job_id: str, feature: Annotated[
    CheckMsaRunningFeature, Depends(check_msa_running_dependency)]) -> CheckJobIsRunningResponse:
    return feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-p2rank-job-running')
async def check_p2rank_job_running(job_id: str, feature: Annotated[
    CheckP2RankRunningFeature, Depends(check_p2rank_running_dependency)]) -> CheckJobIsRunningResponse:
    return feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-umol-job-running')
async def check_umol_job_running(job_id: str, feature: Annotated[
    CheckUmolRunningFeature, Depends(check_umol_running_dependency)]) -> CheckJobIsRunningResponse:
    return feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-result-data-available')
async def check_result_data_available(experiment_id: str,
                                      target_id: str,
                                      ligand_id: str,
                                      job_id: str, feature: Annotated[
            CheckResultDataAvailableFeature, Depends(
                check_result_data_available_dependency)]) -> CheckResultDataAvailableResponse:
    return feature.handle(CheckResultDataAvailableRequest(experiment_id=experiment_id,
                                                          target_id=target_id,
                                                          ligand_id=ligand_id,
                                                          job_id=job_id))


# Docking
@router.post('/run-docking-job')
def perform_docking(request: RunDockingJobRequest, feature: Annotated[
    PredictDockingFeature, Depends(predict_docking_dependency)]) -> RunDockingJobResponse:
    return feature.handle(request)


@router.get('/get-results-list-for-target-ligand')
async def get_results_list_for_target_ligand(experiment_id: str,
                                             target_id: str,
                                             ligand_id: str,
                                             feature: Annotated[GetResultsListForTargetLigandFeature,
                                             Depends(get_results_list_for_target_ligand_feature)]) -> \
        GetResultsListForTargetLigandResponse:
    return feature.handle(GetResultsListForTargetLigandRequest(experiment_id, target_id, ligand_id))


@router.get('/get-all-results-list')
async def get_all_results_list(experiment_id: str,
                                             feature: Annotated[GetAllResultsListFeature,
                                             Depends(
                                                 get_all_results_list_dependency
                                                     )
                                             ]) -> \
        GetAllResultsListResponse:
    return feature.handle(GetAllResultsListRequest(experiment_id=experiment_id))


@router.get('/get-docking-result-data')
async def get_docking_result_data(experiment_id: str,
                                  target_id: str,
                                  ligand_id: str,
                                  job_id: str, feature: Annotated[
            GetDockingResultsFeature, Depends(get_docking_result_dependency)]) -> \
        GetDockingResultDataResponse:
    return feature.handle(GetDockingResultDataRequest(experiment_id, target_id, ligand_id, job_id))
