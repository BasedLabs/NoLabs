import json
from typing import List, Dict, Optional

from fastapi import APIRouter, Depends, UploadFile, File, Form
from typing import Annotated

from nolabs.api_models.experiment import ExperimentMetadataResponse, ChangeExperimentNameRequest, \
    ExperimentMetadataRequest
from nolabs.controllers.drug_discovery.dependencies import (
    get_experiments_feature_dependency,
    delete_experiment_feature_dependency,
    change_experiment_name_dependency,
    upload_target_dependency,
    delete_target_dependency,
    get_targets_list_dependency,
    get_target_data_dependency,
    get_target_ligand_data_dependency,
    upload_target_ligand_dependency,
    delete_target_ligand_dependency,
    get_target_ligands_list_dependency,
    get_binding_pocket_dependency,
    predict_binding_pocket_dependency,
    predict_esmfold_light_dependency,
    predict_esmfold_dependency,
    get_folded_structure_dependency,
    generate_msa_dependency,
    register_docking_job_dependency,
    check_msa_running_dependency,
    check_p2rank_running_dependency,
    check_umol_running_dependency,
    check_result_data_available_dependency,
    get_results_list_for_target_ligand_feature,
    get_all_results_list_dependency,
    predict_umol_docking_dependency,
    get_umol_docking_result_dependency, create_experiment_dependency, get_target_meta_data_dependency,
    get_target_ligand_meta_data_dependency, check_msa_data_available_dependency, check_pocket_data_available_dependency,
    check_folding_data_available_dependency, delete_docking_job_dependency, check_msa_service_health_dependency,
    check_p2rank_service_health_dependency, check_esmfold_service_health_dependency,
    check_esmfold_light_service_health_dependency,
    check_umol_service_health_dependency, get_job_binding_pocket_dependency, set_binding_pocket_dependency,
    get_experiment_metadata_dependency, update_target_name_dependency, get_all_jobs_list_dependency,
    get_jobs_list_for_target_ligand_feature, check_esmfold_running_dependency, check_esmfold_light_running_dependency,
    predict_diffdock_docking_dependency, check_diffdock_running_dependency, check_diffdock_service_health_dependency,
    get_diffdock_ligand_sdf_dependency, get_diffdock_docking_result_dependency, update_diffdock_params_dependency,
    upload_lone_ligand_dependency,
    delete_lone_ligand_dependency, get_lone_ligand_meta_data_dependency, get_lone_ligands_list_dependency,
    get_lone_ligand_data_dependency, get_diffdock_params_dependency, update_docking_params_dependency, predict_rosettafold_dependency,
    check_rosettafold_service_health_dependency, check_rosettafold_running_dependency
)
from nolabs.modules.infrastructure.check_service_health import CheckMsaServiceHealthFeature, \
    CheckP2RankServiceHealthFeature, CheckEsmFoldServiceHealthFeature, \
    CheckUmolServiceHealthFeature, CheckEsmFoldLightServiceHealthFeature, CheckDiffDockServiceHealthFeature, \
    CheckRosettaFoldServiceHealthFeature
from nolabs.modules.drug_discovery.delete_job_feature import DeleteJobFeature
from nolabs.modules.drug_discovery.diffdock_params_management import UpdateDiffDockParamsFeature, \
    GetDiffDockParamsFeature
from nolabs.modules.drug_discovery.docking_params_management import UpdateDockingParamsFeature
from nolabs.modules.drug_discovery.get_experiment_metadata import GetExperimentMetaDataFeature
from nolabs.modules.drug_discovery.lone_ligand_management import UploadLoneLigandFeature, DeleteLoneLigandFeature, \
    GetLoneLigandMetaDataFeature, GetLoneLigandsListFeature, GetLoneLigandDataFeature
from nolabs.modules.drug_discovery.predict_diffdock_docking import PredictDiffDockDockingFeature
from nolabs.modules.drug_discovery.predict_folding import PredictEsmFoldFeature
from nolabs.modules.drug_discovery.predict_rosettafold_folding import PredictRosettaFoldFolding
from nolabs.modules.drug_discovery.result_management import CheckResultDataAvailableFeature, \
    GetAllResultsListFeature, GetResultsListForTargetLigandFeature, CheckMsaDataAvailableFeature, \
    CheckBindingPocketDataAvailableFeature, CheckFoldingDataAvailableFeature, GetBJobBindingPocketDataFeature, \
    GetAllJobsListFeature, GetJobsListForTargetLigandFeature
from nolabs.modules.drug_discovery.set_binding_pocket import SetBindingPocketFeature
from nolabs.modules.drug_discovery.target_management import UploadTargetFeature, DeleteTargetFeature, \
    GetTargetsListFeature, GetTargetDataFeature, GetTargetMetaDataFeature, UpdateTargetNameFeature
from nolabs.modules.drug_discovery.get_binding_pocket import GetBindingPocketFeature
from nolabs.modules.drug_discovery.predict_binding_pocket import PredictBindingPocketFeature
from nolabs.modules.drug_discovery.predict_light_folding import PredictEsmFoldLightFeature
from nolabs.modules.drug_discovery.get_folding import GetFoldedStructureFeature
from nolabs.modules.drug_discovery.generate_msa import GenerateMsaFeature
from nolabs.modules.drug_discovery.target_ligand_management import UploadTargetLigandFeature, DeleteTargetLigandFeature, \
    GetTargetLigandsListFeature, GetTargetLigandDataFeature, GetTargetLigandMetaDataFeature
from nolabs.modules.drug_discovery.register_docking_job import RegisterDockingJobFeature
from nolabs.modules.drug_discovery.predict_umol_docking import PredictUmolDockingFeature
from nolabs.modules.experiment.create_experiment import CreateExperimentFeature
from nolabs.modules.drug_discovery.get_umol_results import GetUmolDockingResultsFeature
from nolabs.modules.drug_discovery.get_diffdock_results import GetDiffDockLigandSdfFeature, \
    GetDiffDockDockingResultsFeature
from nolabs.modules.drug_discovery.progress_management import CheckMsaRunningFeature, \
    CheckP2RankRunningFeature, CheckUmolRunningFeature, CheckEsmFoldRunningFeature, CheckEsmFoldLightRunningFeature, \
    CheckDiffDockRunningFeature, CheckRosettaFoldRunningFeature
from nolabs.modules.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.modules.experiment.change_experiment_name import ChangeExperimentNameFeature

from nolabs.api_models.drug_discovery import (
    UploadTargetRequest,
    UploadTargetResponse,
    DeleteTargetRequest,
    DeleteTargetResponse,
    GetTargetsListRequest,
    GetTargetDataRequest,
    GetTargetDataResponse,
    UploadTargetLigandRequest,
    UploadTargetLigandResponse,
    DeleteTargetLigandRequest,
    DeleteTargetLigandResponse,
    GetTargetLigandsListRequest,
    GetTargetLigandDataResponse,
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
    RunUmolDockingJobRequest,
    RunUmolDockingJobResponse,
    CheckResultDataAvailableRequest,
    CheckResultDataAvailableResponse,
    GetDockingResultDataRequest,
    GetUmolDockingResultDataResponse,
    TargetMetaData,
    LigandMetaData,
    GetTargetLigandDataRequest,
    GetResultsListForTargetLigandRequest,
    GetResultsListForTargetLigandResponse,
    GetAllResultsListRequest,
    GetAllResultsListResponse, GetTargetMetaDataResponse, GetTargetMetaDataRequest, GetTargetLigandMetaDataResponse,
    GetTargetLigandMetaDataRequest, CheckMsaDataAvailableResponse, CheckMsaDataAvailableRequest,
    CheckPocketDataAvailableResponse, CheckPocketDataAvailableRequest, CheckFoldingDataAvailableResponse,
    CheckFoldingDataAvailableRequest, DeleteDockingJobRequest, DeleteDockingJobResponse, CheckServiceHealthyResponse,
    GetJobBindingPocketDataRequest, GetJobBindingPocketDataResponse, SetTargetBindingPocketRequest,
    UpdateTargetNameRequest, GetAllJobsListResponse, GetAllJobsListRequest, GetJobsListForTargetLigandResponse,
    GetJobsListForTargetLigandRequest, RunDiffDockDockingJobResponse, RunDiffDockDockingJobRequest,
    GetDiffDockLigandSdfResponse, GetDiffDockLigandSdfRequest, GetDiffDockDockingResultDataResponse,
    UpdateDiffDockParamsRequest, GetDiffDockParamsResponse, GetDiffDockParamsRequest, GetDockingParamsResponse,
    UpdateDockingParamsRequest, UploadLoneLigandResponse, UploadLoneLigandRequest, DeleteLoneLigandResponse,
    DeleteLoneLigandRequest, GetLoneLigandMetaDataResponse, GetLoneLigandMetaDataRequest, GetLoneLigandsListRequest,
    GetLoneLigandDataResponse, GetLoneLigandDataRequest
)
from nolabs.modules.experiment.get_experiments import GetExperimentsFeature

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


@router.get('/get-experiment-metadata')
async def get_experiment_metadata(experiment_id: str,
                                  feature: Annotated[
                                      GetExperimentMetaDataFeature, Depends(
                                          get_experiment_metadata_dependency)]) -> ExperimentMetadataResponse:
    return feature.handle(ExperimentMetadataRequest(experiment_id))


@router.post('/change-experiment-name')
async def change_experiment_name(experiment_id: str, experiment_name: str, feature: Annotated[
    ChangeExperimentNameFeature, Depends(change_experiment_name_dependency)]) -> ExperimentMetadataResponse:
    return feature.handle(ChangeExperimentNameRequest(experiment_id, experiment_name))


# Target Management
@router.post('/upload-target')
async def upload_target(feature: Annotated[UploadTargetFeature, Depends(upload_target_dependency)],
                        experiment_id: str = Form(),
                        fasta: UploadFile = File(),
                        metadata: Optional[str] = Form(default=None),  # Now expects a JSON string for metadata
                        ) -> UploadTargetResponse:
    metadata_dict = json.loads(metadata) if metadata else None
    return feature.handle(UploadTargetRequest(experiment_id, fasta, metadata_dict))


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


@router.get('/get-target-meta-data')
async def get_target_meta_data(feature: Annotated[
    GetTargetMetaDataFeature, Depends(get_target_meta_data_dependency)],
                               experiment_id: str,
                               target_id: str,
                               ) -> GetTargetMetaDataResponse:
    result = feature.handle(GetTargetMetaDataRequest(experiment_id, target_id))
    return GetTargetMetaDataResponse(target_id=result.target_id, target_name=result.target_name)


@router.post('/update-target-name')
async def update_target_name(feature: Annotated[
    UpdateTargetNameFeature, Depends(update_target_name_dependency)],
                             experiment_id: str,
                             target_id: str,
                             target_name: str
                             ) -> GetTargetMetaDataResponse:
    result = feature.handle(UpdateTargetNameRequest(experiment_id, target_id, target_name))
    return GetTargetMetaDataResponse(target_id=result.target_id, target_name=result.target_name)


@router.get('/get-target-data')
async def get_target_data(feature: Annotated[
    GetTargetDataFeature, Depends(get_target_data_dependency)],
                          experiment_id: str,
                          target_id: str,
                          ) -> GetTargetDataResponse:
    result = feature.handle(GetTargetDataRequest(experiment_id, target_id))
    return GetTargetDataResponse(protein_sequence=result.protein_sequence, protein_pdb=result.protein_pdb)


# Ligand Management
@router.post('/upload-ligand-to-target')
async def upload_ligand_to_target(feature: Annotated[UploadTargetLigandFeature, Depends(upload_target_ligand_dependency)],
                        experiment_id: str = Form(),
                        target_id: str = Form(),
                        sdf_file: UploadFile = File()
                        ) -> UploadTargetLigandResponse:
    return feature.handle(UploadTargetLigandRequest(experiment_id, target_id, sdf_file))

@router.post('/upload-ligand-to-experiment')
async def upload_ligand_to_experiment(feature: Annotated[UploadLoneLigandFeature, Depends(upload_lone_ligand_dependency)],
                        experiment_id: str = Form(),
                        sdf_file: UploadFile = File(),
                        metadata: Optional[str] = Form(default=None)
                        ) -> UploadLoneLigandResponse:
    metadata_dict = json.loads(metadata) if metadata else None
    return feature.handle(UploadLoneLigandRequest(experiment_id, sdf_file, metadata_dict))

@router.delete('/delete-ligand-from-target')
async def delete_ligand_from_target(feature: Annotated[
    DeleteTargetLigandFeature, Depends(delete_target_ligand_dependency)],
                        experiment_id: str,
                        target_id: str,
                        ligand_id: str) -> DeleteTargetLigandResponse:
    return feature.handle(DeleteTargetLigandRequest(experiment_id, target_id, ligand_id))

@router.delete('/delete-ligand-from-experiment')
async def delete_ligand_from_experiment(feature: Annotated[
    DeleteLoneLigandFeature, Depends(delete_lone_ligand_dependency)],
                        experiment_id: str,
                        ligand_id: str) -> DeleteLoneLigandResponse:
    return feature.handle(DeleteLoneLigandRequest(experiment_id, ligand_id))

@router.get('/get-ligand-meta-data-for-target')
async def get_ligand_meta_data_for_target(feature: Annotated[
    GetTargetLigandMetaDataFeature, Depends(get_target_ligand_meta_data_dependency)],
                          experiment_id: str,
                          target_id: str,
                          ligand_id: str) -> GetTargetLigandMetaDataResponse:
    return feature.handle(GetTargetLigandMetaDataRequest(experiment_id, target_id, ligand_id))

@router.get('/get-lone-ligand-meta-data')
async def get_lone_ligand_meta_data(feature: Annotated[
    GetLoneLigandMetaDataFeature, Depends(get_lone_ligand_meta_data_dependency)],
                          experiment_id: str,
                          ligand_id: str) -> GetLoneLigandMetaDataResponse:
    return feature.handle(GetLoneLigandMetaDataRequest(experiment_id, ligand_id))


@router.get('/get-ligand-data-for-target')
async def get_ligand_data_for_target(feature: Annotated[
    GetTargetLigandDataFeature, Depends(get_target_ligand_data_dependency)],
                          experiment_id: str,
                          target_id: str,
                          ligand_id: str) -> GetTargetLigandDataResponse:
    return feature.handle(GetTargetLigandDataRequest(experiment_id, target_id, ligand_id))

@router.get('/get-lone-ligand-data')
async def get_lone_ligand_data(feature: Annotated[
    GetLoneLigandDataFeature, Depends(get_lone_ligand_data_dependency)],
                          experiment_id: str,
                          ligand_id: str) -> GetLoneLigandDataResponse:
    return feature.handle(GetLoneLigandDataRequest(experiment_id, ligand_id))



@router.get('/get-ligands-list-for-target')
async def get_ligands_list_for_target(feature: Annotated[
    GetTargetLigandsListFeature, Depends(get_target_ligands_list_dependency)],
                           experiment_id: str,
                           target_id: str) -> List[LigandMetaData]:
    return feature.handle(GetTargetLigandsListRequest(experiment_id, target_id)).ligands

@router.get('/get-lone-ligands-list')
async def get_lone_ligands_list(feature: Annotated[
    GetLoneLigandsListFeature, Depends(get_lone_ligands_list_dependency)],
                           experiment_id: str) -> List[LigandMetaData]:
    return feature.handle(GetLoneLigandsListRequest(experiment_id)).ligands

# Binding Pocket Management
@router.get('/get-target-binding-pocket')
async def get_target_binding_pocket(feature: Annotated[GetBindingPocketFeature, Depends(get_binding_pocket_dependency)],
                                    experiment_id: str,
                                    target_id: str,
                                    ) -> GetTargetBindingPocketResponse:
    return feature.handle(GetTargetBindingPocketRequest(experiment_id, target_id))


@router.post('/set-target-binding-pocket')
async def set_target_binding_pocket(feature: Annotated[SetBindingPocketFeature, Depends(set_binding_pocket_dependency)],
                                    experiment_id: str,
                                    target_id: str,
                                    pocket_ids: List[int]
                                    ):
    return feature.handle(SetTargetBindingPocketRequest(experiment_id, target_id, pocket_ids))


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
                               folding_method: str,
                               feature: Annotated[GetFoldedStructureFeature,
                               Depends(get_folded_structure_dependency)]) -> GetFoldingResponse:
    return feature.handle(GetFoldingRequest(experiment_id, target_id, folding_method))


@router.post('/predict-esmfold-light')
async def predict_esmfold_light(feature: Annotated[
    PredictEsmFoldLightFeature, Depends(predict_esmfold_light_dependency)],
                                experiment_id: str,
                                target_id: str
                                ) -> PredictFoldingResponse:
    return feature.handle(PredictFoldingRequest(experiment_id, target_id))


@router.post('/predict-rosettafold')
async def predict_rosettafold(feature: Annotated[
    PredictRosettaFoldFolding, Depends(predict_rosettafold_dependency)],
                              experiment_id: str,
                              target_id: str
                              ) -> PredictFoldingResponse:
    return feature.handle(PredictFoldingRequest(experiment_id, target_id))


@router.post('/predict-esmfold')
async def predict_esmfold(feature: Annotated[
    PredictEsmFoldFeature, Depends(predict_esmfold_dependency)],
                          experiment_id: str,
                          target_id: str
                          ) -> PredictFoldingResponse:
    return feature.handle(PredictFoldingRequest(experiment_id, target_id))


@router.post('/register-docking-job')
async def register_docking_job(experiment_id: str,
                               target_id: str,
                               ligand_id: str,
                               folding_method: str,
                               feature: Annotated[RegisterDockingJobFeature
                               , Depends(register_docking_job_dependency)]) -> RegisterDockingJobResponse:
    return feature.handle(RegisterDockingJobRequest(experiment_id, target_id, ligand_id, folding_method))


@router.delete('/delete-docking-job')
async def delete_docking_job(experiment_id: str,
                             target_id: str,
                             ligand_id: str,
                             job_id: str,
                             feature: Annotated[DeleteJobFeature
                             , Depends(delete_docking_job_dependency)]) -> DeleteDockingJobResponse:
    return feature.handle(DeleteDockingJobRequest(experiment_id, target_id, ligand_id, job_id))


@router.get('/check-msa-data-available')
async def check_msa_data_available(experiment_id: str,
                                   target_id: str, feature: Annotated[
            CheckMsaDataAvailableFeature, Depends(
                check_msa_data_available_dependency)]) -> CheckMsaDataAvailableResponse:
    return feature.handle(CheckMsaDataAvailableRequest(experiment_id, target_id))


@router.get('/check-msa-service-health')
async def check_msa_service_health(feature: Annotated[
    CheckMsaServiceHealthFeature, Depends(check_msa_service_health_dependency)]) -> CheckServiceHealthyResponse:
    return await feature.handle()


@router.get('/check-msa-job-running')
async def check_msa_job_running(job_id: str, feature: Annotated[
    CheckMsaRunningFeature, Depends(check_msa_running_dependency)]) -> CheckJobIsRunningResponse:
    return await feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-rosettafold-job-running')
async def check_rosettafold_job_running(job_id: str, feature: Annotated[
    CheckRosettaFoldRunningFeature, Depends(check_rosettafold_running_dependency)]) -> CheckJobIsRunningResponse:
    return await feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-pocket-data-available')
async def check_pocket_data_available(experiment_id: str,
                                      target_id: str,
                                      ligand_id: str,
                                      job_id: str,
                                      feature: Annotated[
                                          CheckBindingPocketDataAvailableFeature, Depends(
                                              check_pocket_data_available_dependency)]) -> CheckPocketDataAvailableResponse:
    return feature.handle(CheckPocketDataAvailableRequest(experiment_id, target_id, ligand_id, job_id))


@router.get('/get-job-pocket-data')
async def get_job_pocket_data(experiment_id: str,
                              target_id: str,
                              ligand_id: str,
                              job_id: str,
                              feature: Annotated[
                                  GetBJobBindingPocketDataFeature, Depends(
                                      get_job_binding_pocket_dependency)]) -> GetJobBindingPocketDataResponse:
    return feature.handle(GetJobBindingPocketDataRequest(experiment_id, target_id, ligand_id, job_id))


@router.get('/check-p2rank-service-health')
async def check_p2rank_service_health(feature: Annotated[
    CheckP2RankServiceHealthFeature, Depends(check_p2rank_service_health_dependency)]) -> CheckServiceHealthyResponse:
    return await feature.handle()


@router.get('/check-p2rank-job-running')
async def check_p2rank_job_running(job_id: str, feature: Annotated[
    CheckP2RankRunningFeature, Depends(check_p2rank_running_dependency)]) -> CheckJobIsRunningResponse:
    return await feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-folding-data-available')
async def check_folding_data_available(experiment_id: str,
                                       target_id: str,
                                       folding_method: str,
                                       feature: Annotated[
                                           CheckFoldingDataAvailableFeature, Depends(
                                               check_folding_data_available_dependency)]) -> CheckFoldingDataAvailableResponse:
    return feature.handle(CheckFoldingDataAvailableRequest(experiment_id, target_id, folding_method))


@router.get('/check-esmfold-service-health')
async def check_esmfold_service_health(feature: Annotated[
    CheckEsmFoldServiceHealthFeature, Depends(check_esmfold_service_health_dependency)]) -> CheckServiceHealthyResponse:
    return await feature.handle()


@router.get('/check-esmfold-light-service-health')
async def check_esmfold_light_service_health(feature: Annotated[
    CheckEsmFoldLightServiceHealthFeature, Depends(
        check_esmfold_light_service_health_dependency)]) -> CheckServiceHealthyResponse:
    return await feature.handle()


@router.get('/check-esmfold-job-running')
async def check_esmfold_job_running(job_id: str, feature: Annotated[
    CheckEsmFoldRunningFeature, Depends(check_esmfold_running_dependency)]) -> CheckJobIsRunningResponse:
    return await feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-esmfold-light-job-running')
async def check_esmfold_light_job_running(job_id: str, feature: Annotated[
    CheckEsmFoldLightRunningFeature, Depends(check_esmfold_light_running_dependency)]) -> CheckJobIsRunningResponse:
    return await feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-umol-service-health')
async def check_umol_service_health(feature: Annotated[
    CheckUmolServiceHealthFeature, Depends(check_umol_service_health_dependency)]) -> CheckServiceHealthyResponse:
    return await feature.handle()


@router.get('/check-rosettafold-service-health')
async def rosettafold_service_health(feature: Annotated[
    CheckRosettaFoldServiceHealthFeature, Depends(
        check_rosettafold_service_health_dependency)]) -> CheckServiceHealthyResponse:
    return await feature.handle()


@router.get('/check-umol-job-running')
async def check_umol_job_running(job_id: str, feature: Annotated[
    CheckUmolRunningFeature, Depends(check_umol_running_dependency)]) -> CheckJobIsRunningResponse:
    return await feature.handle(CheckJobIsRunningRequest(job_id=job_id))


@router.get('/check-diffdock-service-health')
async def check_diffdock_service_health(feature: Annotated[
    CheckDiffDockServiceHealthFeature, Depends(
        check_diffdock_service_health_dependency)]) -> CheckServiceHealthyResponse:
    return await feature.handle()


@router.get('/check-diffdock-job-running')
async def check_diffdock_job_running(job_id: str, feature: Annotated[
    CheckDiffDockRunningFeature, Depends(check_diffdock_running_dependency)]) -> CheckJobIsRunningResponse:
    return await feature.handle(CheckJobIsRunningRequest(job_id=job_id))


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
@router.post('/run-umol-docking-job')
def run_umol_docking(experiment_id: str,
                     target_id: str,
                     ligand_id: str,
                     job_id: str
                     , feature: Annotated[
            PredictUmolDockingFeature, Depends(predict_umol_docking_dependency)]) -> RunUmolDockingJobResponse:
    return feature.handle(RunUmolDockingJobRequest(experiment_id=experiment_id,
                                                   target_id=target_id,
                                                   ligand_id=ligand_id, job_id=job_id))


@router.post('/run-diffdock-docking-job')
def run_diffdock_docking(experiment_id: str,
                         target_id: str,
                         ligand_id: str,
                         job_id: str
                         , feature: Annotated[
            PredictDiffDockDockingFeature, Depends(
                predict_diffdock_docking_dependency)]) -> RunDiffDockDockingJobResponse:
    return feature.handle(RunDiffDockDockingJobRequest(experiment_id=experiment_id,
                                                       target_id=target_id,
                                                       ligand_id=ligand_id,
                                                       job_id=job_id))


@router.get('/get-diffdock-params')
def get_diffdock_params(experiment_id: str,
                        target_id: str,
                        ligand_id: str,
                        job_id: str,
                        feature: Annotated[
                            GetDiffDockParamsFeature, Depends(
                                get_diffdock_params_dependency)]) -> GetDiffDockParamsResponse:
    return feature.handle(GetDiffDockParamsRequest(experiment_id=experiment_id,
                                                   target_id=target_id,
                                                   ligand_id=ligand_id,
                                                   job_id=job_id))


@router.post('/update-docking-params')
def update_docking_params(experiment_id: str,
                          target_id: str,
                          ligand_id: str,
                          job_id: str,
                          folding_method: str,
                          docking_method: str,
                          feature: Annotated[
                              UpdateDockingParamsFeature, Depends(
                                  update_docking_params_dependency)]) -> GetDockingParamsResponse:
    return feature.handle(UpdateDockingParamsRequest(experiment_id=experiment_id,
                                                     target_id=target_id,
                                                     ligand_id=ligand_id,
                                                     job_id=job_id,
                                                     folfing_method=folding_method,
                                                     docking_method=docking_method))


@router.post('/update-diffdock-params')
def update_diffdock_params(experiment_id: str,
                           target_id: str,
                           ligand_id: str,
                           job_id: str,
                           samples_per_complex: int,
                           feature: Annotated[
                               UpdateDiffDockParamsFeature, Depends(
                                   update_diffdock_params_dependency)]) -> GetDiffDockParamsResponse:
    return feature.handle(UpdateDiffDockParamsRequest(experiment_id=experiment_id,
                                                      target_id=target_id,
                                                      ligand_id=ligand_id,
                                                      job_id=job_id,
                                                      samples_per_complex=samples_per_complex))


@router.get('/get-results-list-for-target-ligand')
async def get_results_list_for_target_ligand(experiment_id: str,
                                             target_id: str,
                                             ligand_id: str,
                                             feature: Annotated[GetResultsListForTargetLigandFeature,
                                             Depends(get_results_list_for_target_ligand_feature)]) -> \
        GetResultsListForTargetLigandResponse:
    return feature.handle(GetResultsListForTargetLigandRequest(experiment_id, target_id, ligand_id))


@router.get('/get-jobs-list-for-target-ligand')
async def get_jobs_list_for_target_ligand(experiment_id: str,
                                          target_id: str,
                                          ligand_id: str,
                                          feature: Annotated[GetJobsListForTargetLigandFeature,
                                          Depends(get_jobs_list_for_target_ligand_feature)]) -> \
        GetJobsListForTargetLigandResponse:
    return feature.handle(GetJobsListForTargetLigandRequest(experiment_id, target_id, ligand_id))


@router.get('/get-all-results-list')
async def get_all_results_list(experiment_id: str,
                               feature: Annotated[GetAllResultsListFeature,
                               Depends(
                                   get_all_results_list_dependency
                               )
                               ]) -> \
        GetAllResultsListResponse:
    return feature.handle(GetAllResultsListRequest(experiment_id=experiment_id))


@router.get('/get-all-jobs-list')
async def get_all_jobs_list(experiment_id: str,
                            feature: Annotated[GetAllJobsListFeature,
                            Depends(
                                get_all_jobs_list_dependency
                            )
                            ]) -> \
        GetAllJobsListResponse:
    return feature.handle(GetAllJobsListRequest(experiment_id=experiment_id))


@router.get('/get-umol-docking-result-data')
async def get_umol_docking_result_data(experiment_id: str,
                                       target_id: str,
                                       ligand_id: str,
                                       job_id: str, feature: Annotated[
            GetUmolDockingResultsFeature, Depends(get_umol_docking_result_dependency)]) -> \
        GetUmolDockingResultDataResponse:
    return feature.handle(GetDockingResultDataRequest(experiment_id, target_id, ligand_id, job_id))


@router.get('/get-diffdock-docking-result-data')
async def get_diffdock_docking_result_data(experiment_id: str,
                                           target_id: str,
                                           ligand_id: str,
                                           job_id: str, feature: Annotated[
            GetDiffDockDockingResultsFeature, Depends(get_diffdock_docking_result_dependency)]) -> \
        GetDiffDockDockingResultDataResponse:
    return feature.handle(GetDockingResultDataRequest(experiment_id, target_id, ligand_id, job_id))


@router.get('/get-diffdock-ligand_sdf')
async def get_diffdock_ligand_sdf(experiment_id: str,
                                  target_id: str,
                                  ligand_id: str,
                                  job_id: str,
                                  sdf_file_name: str,
                                  feature: Annotated[
                                      GetDiffDockLigandSdfFeature, Depends(get_diffdock_ligand_sdf_dependency)]) -> \
        GetDiffDockLigandSdfResponse:
    return feature.handle(GetDiffDockLigandSdfRequest(experiment_id, target_id, ligand_id, job_id, sdf_file_name))
