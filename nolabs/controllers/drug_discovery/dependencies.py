from typing import Annotated

from fastapi import Depends

from nolabs.controllers.common_dependencies import settings_dependency
from nolabs.modules.infrastructure.check_service_health import CheckMsaServiceHealthFeature, \
    CheckP2RankServiceHealthFeature, CheckUmolServiceHealthFeature, \
    CheckEsmFoldServiceHealthFeature, CheckEsmFoldLightServiceHealthFeature, CheckDiffDockServiceHealthFeature, \
    CheckRosettaFoldServiceHealthFeature
from nolabs.modules.drug_discovery.delete_job_feature import DeleteJobFeature
from nolabs.modules.drug_discovery.diffdock_params_management import GetDiffDockParamsFeature, \
    UpdateDiffDockParamsFeature
from nolabs.modules.drug_discovery.docking_params_management import UpdateDockingParamsFeature
from nolabs.modules.drug_discovery.get_experiment_metadata import GetExperimentMetaDataFeature
from nolabs.modules.drug_discovery.lone_ligand_management import UploadLoneLigandFeature, DeleteLoneLigandFeature, \
    GetLoneLigandsListFeature, GetLoneLigandMetaDataFeature, GetLoneLigandDataFeature
from nolabs.modules.drug_discovery.predict_diffdock_docking import PredictDiffDockDockingFeature
from nolabs.modules.drug_discovery.predict_folding import PredictEsmFoldFeature
from nolabs.modules.drug_discovery.predict_rosettafold_folding import PredictRosettaFoldFolding
from nolabs.modules.drug_discovery.result_management import CheckResultDataAvailableFeature, \
    GetAllResultsListFeature, GetResultsListForTargetLigandFeature, CheckMsaDataAvailableFeature, \
    CheckBindingPocketDataAvailableFeature, CheckFoldingDataAvailableFeature, GetBJobBindingPocketDataFeature, \
    GetAllJobsListFeature, GetJobsListForTargetLigandFeature
from nolabs.modules.drug_discovery.set_binding_pocket import SetBindingPocketFeature
from nolabs.modules.experiment.create_experiment import CreateExperimentFeature
from nolabs.modules.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.modules.experiment.change_experiment_name import ChangeExperimentNameFeature
from nolabs.modules.drug_discovery.generate_msa import GenerateMsaFeature
from nolabs.modules.drug_discovery.services.file_management import FileManagement
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement
from nolabs.modules.drug_discovery.target_management import UploadTargetFeature, DeleteTargetFeature, \
    GetTargetsListFeature, GetTargetDataFeature, GetTargetMetaDataFeature, UpdateTargetNameFeature
from nolabs.modules.drug_discovery.target_ligand_management import UploadTargetLigandFeature, \
    DeleteTargetLigandFeature, \
    GetTargetLigandsListFeature, GetTargetLigandDataFeature, GetTargetLigandMetaDataFeature
from nolabs.modules.drug_discovery.get_binding_pocket import GetBindingPocketFeature
from nolabs.modules.drug_discovery.predict_binding_pocket import PredictBindingPocketFeature
from nolabs.modules.drug_discovery.predict_light_folding import PredictEsmFoldLightFeature
from nolabs.modules.drug_discovery.get_folding import GetFoldedStructureFeature
from nolabs.modules.drug_discovery.register_docking_job import RegisterDockingJobFeature
from nolabs.modules.drug_discovery.progress_management import CheckMsaRunningFeature, \
    CheckP2RankRunningFeature, CheckUmolRunningFeature, CheckEsmFoldRunningFeature, CheckEsmFoldLightRunningFeature, \
    CheckDiffDockRunningFeature, CheckRosettaFoldRunningFeature
from nolabs.modules.drug_discovery.predict_umol_docking import PredictUmolDockingFeature
from nolabs.modules.drug_discovery.get_umol_results import GetUmolDockingResultsFeature
from nolabs.modules.drug_discovery.get_diffdock_results import GetDiffDockDockingResultsFeature, \
    GetDiffDockLigandSdfFeature
from nolabs.modules.experiment.get_experiments import GetExperimentsFeature
from nolabs.infrastructure.settings import Settings


def file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> FileManagement:
    return FileManagement(settings=settings)


def get_experiments_feature_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> GetExperimentsFeature:
    return GetExperimentsFeature(file_management=file_management)


def get_experiment_metadata_dependency(
        file_management: Annotated[
            FileManagement, Depends(file_management_dependency)]) -> GetExperimentMetaDataFeature:
    return GetExperimentMetaDataFeature(file_management=file_management)


def delete_experiment_feature_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> DeleteExperimentFeature:
    return DeleteExperimentFeature(file_management=file_management)


def change_experiment_name_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]
) -> ChangeExperimentNameFeature:
    return ChangeExperimentNameFeature(
        file_management=file_management
    )


def target_file_management_dependency(
        settings: Annotated[Settings, Depends(settings_dependency)]) -> TargetsFileManagement:
    return TargetsFileManagement(settings=settings)


def ligand_file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                                      experiments_file_management: Annotated[
                                          FileManagement, Depends(file_management_dependency)],
                                      target_file_management: Annotated[TargetsFileManagement,
                                      Depends(target_file_management_dependency)]) -> LigandsFileManagement:
    return LigandsFileManagement(settings=settings, experiments_file_management=experiments_file_management,
                                 targets_file_management=target_file_management)


def result_file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                                      ligand_file_management: Annotated[LigandsFileManagement,
                                      Depends(ligand_file_management_dependency)]) -> ResultsFileManagement:
    return ResultsFileManagement(settings=settings, ligand_file_management=ligand_file_management)


def generate_msa_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                            target_file_management: Annotated[TargetsFileManagement,
                            Depends(target_file_management_dependency)]) -> GenerateMsaFeature:
    return GenerateMsaFeature(file_management=target_file_management, settings=settings)


def upload_target_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> UploadTargetFeature:
    return UploadTargetFeature(file_management=target_file_management)


def delete_target_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> DeleteTargetFeature:
    return DeleteTargetFeature(file_management=target_file_management)


def get_targets_list_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> GetTargetsListFeature:
    return GetTargetsListFeature(file_management=target_file_management)


def get_target_meta_data_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> GetTargetMetaDataFeature:
    return GetTargetMetaDataFeature(file_management=target_file_management)


def update_target_name_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> UpdateTargetNameFeature:
    return UpdateTargetNameFeature(file_management=target_file_management)


def get_target_data_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> GetTargetDataFeature:
    return GetTargetDataFeature(file_management=target_file_management)


def get_binding_pocket_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> GetBindingPocketFeature:
    return GetBindingPocketFeature(file_management=target_file_management)


def set_binding_pocket_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> SetBindingPocketFeature:
    return SetBindingPocketFeature(file_management=target_file_management)


def predict_binding_pocket_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                                      target_file_management: Annotated[TargetsFileManagement,
                                      Depends(target_file_management_dependency)]) -> PredictBindingPocketFeature:
    return PredictBindingPocketFeature(file_management=target_file_management, settings=settings)


def get_folded_structure_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> GetFoldedStructureFeature:
    return GetFoldedStructureFeature(file_management=target_file_management)


def predict_esmfold_light_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)],
                                     settings: Annotated[Settings,
                                     Depends(settings_dependency)]) -> PredictEsmFoldLightFeature:
    return PredictEsmFoldLightFeature(file_management=target_file_management, settings=settings)


def predict_esmfold_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)],
                               settings: Annotated[Settings,
                               Depends(settings_dependency)]) -> PredictEsmFoldFeature:
    return PredictEsmFoldFeature(file_management=target_file_management, settings=settings)

def predict_rosettafold_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)],
                               settings: Annotated[Settings,
                               Depends(settings_dependency)]) -> PredictRosettaFoldFolding:
    return PredictRosettaFoldFolding(file_management=target_file_management, settings=settings)


def upload_target_ligand_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> UploadTargetLigandFeature:
    return UploadTargetLigandFeature(file_management=ligand_file_management)


def upload_lone_ligand_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> UploadLoneLigandFeature:
    return UploadLoneLigandFeature(file_management=ligand_file_management)


def delete_target_ligand_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> DeleteTargetLigandFeature:
    return DeleteTargetLigandFeature(file_management=ligand_file_management)

def delete_lone_ligand_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> DeleteLoneLigandFeature:
    return DeleteLoneLigandFeature(file_management=ligand_file_management)

def get_target_ligands_list_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> GetTargetLigandsListFeature:
    return GetTargetLigandsListFeature(file_management=ligand_file_management)

def get_lone_ligands_list_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> GetLoneLigandsListFeature:
    return GetLoneLigandsListFeature(file_management=ligand_file_management)


def get_target_ligand_meta_data_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> GetTargetLigandMetaDataFeature:
    return GetTargetLigandMetaDataFeature(file_management=ligand_file_management)

def get_lone_ligand_meta_data_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> GetLoneLigandMetaDataFeature:
    return GetLoneLigandMetaDataFeature(file_management=ligand_file_management)


def get_target_ligand_data_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> GetTargetLigandDataFeature:
    return GetTargetLigandDataFeature(file_management=ligand_file_management)


def get_lone_ligand_data_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> GetLoneLigandDataFeature:
    return GetLoneLigandDataFeature(file_management=ligand_file_management)

def register_docking_job_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]):
    return RegisterDockingJobFeature(file_management=result_file_management)

def delete_docking_job_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]):
    return DeleteJobFeature(file_management=result_file_management)


def check_msa_data_available_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> CheckMsaDataAvailableFeature:
    return CheckMsaDataAvailableFeature(file_management=target_file_management)


def check_msa_service_health_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckMsaServiceHealthFeature:
    return CheckMsaServiceHealthFeature(settings=settings)


def check_msa_running_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckMsaRunningFeature:
    return CheckMsaRunningFeature(settings=settings)


def check_rosettafold_running_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckRosettaFoldRunningFeature:
    return CheckRosettaFoldRunningFeature(settings=settings)


def check_pocket_data_available_dependency(results_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> CheckBindingPocketDataAvailableFeature:
    return CheckBindingPocketDataAvailableFeature(file_management=results_file_management)


def get_job_binding_pocket_dependency(results_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> GetBJobBindingPocketDataFeature:
    return GetBJobBindingPocketDataFeature(file_management=results_file_management)


def check_p2rank_service_health_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckP2RankServiceHealthFeature:
    return CheckP2RankServiceHealthFeature(settings=settings)


def check_p2rank_running_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckP2RankRunningFeature:
    return CheckP2RankRunningFeature(settings=settings)


def check_folding_data_available_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> CheckFoldingDataAvailableFeature:
    return CheckFoldingDataAvailableFeature(file_management=target_file_management)


def check_esmfold_service_health_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckEsmFoldServiceHealthFeature:
    return CheckEsmFoldServiceHealthFeature(settings=settings)


def check_esmfold_light_service_health_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckEsmFoldLightServiceHealthFeature:
    return CheckEsmFoldLightServiceHealthFeature(settings=settings)


def check_esmfold_running_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckEsmFoldRunningFeature:
    return CheckEsmFoldRunningFeature(settings=settings)


def check_esmfold_light_running_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckEsmFoldLightRunningFeature:
    return CheckEsmFoldLightRunningFeature(settings=settings)


def check_umol_service_health_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckUmolServiceHealthFeature:
    return CheckUmolServiceHealthFeature(settings=settings)


def check_rosettafold_service_health_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckRosettaFoldServiceHealthFeature:
    return CheckRosettaFoldServiceHealthFeature(settings=settings)


def check_umol_running_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckUmolRunningFeature:
    return CheckUmolRunningFeature(settings=settings)


def check_diffdock_service_health_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckDiffDockServiceHealthFeature:
    return CheckDiffDockServiceHealthFeature(settings=settings)


def check_diffdock_running_dependency(settings: Annotated[Settings,
Depends(settings_dependency)]) -> CheckDiffDockRunningFeature:
    return CheckDiffDockRunningFeature(settings=settings)


def check_result_data_available_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> CheckResultDataAvailableFeature:
    return CheckResultDataAvailableFeature(file_management=result_file_management)


def get_all_results_list_dependency(
        target_file_management: Annotated[TargetsFileManagement,
        Depends(target_file_management_dependency)], ligand_file_management: Annotated[LigandsFileManagement,
        Depends(ligand_file_management_dependency)], result_file_management: Annotated[ResultsFileManagement,
        Depends(result_file_management_dependency)]) -> GetAllResultsListFeature:
    return GetAllResultsListFeature(target_file_management=target_file_management,
                                    ligand_file_management=ligand_file_management,
                                    results_file_management=result_file_management)


def get_all_jobs_list_dependency(
        target_file_management: Annotated[TargetsFileManagement,
        Depends(target_file_management_dependency)], ligand_file_management: Annotated[LigandsFileManagement,
        Depends(ligand_file_management_dependency)], result_file_management: Annotated[ResultsFileManagement,
        Depends(result_file_management_dependency)]) -> GetAllJobsListFeature:
    return GetAllJobsListFeature(target_file_management=target_file_management,
                                 ligand_file_management=ligand_file_management,
                                 results_file_management=result_file_management)


def get_results_list_for_target_ligand_feature(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> GetResultsListForTargetLigandFeature:
    return GetResultsListForTargetLigandFeature(results_file_management=result_file_management)


def get_jobs_list_for_target_ligand_feature(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> GetJobsListForTargetLigandFeature:
    return GetJobsListForTargetLigandFeature(results_file_management=result_file_management)


def predict_umol_docking_dependency(
        settings: Annotated[Settings, Depends(settings_dependency)],
        target_file_management: Annotated[TargetsFileManagement,
        Depends(target_file_management_dependency)], ligand_file_management: Annotated[LigandsFileManagement,
        Depends(ligand_file_management_dependency)], result_file_management: Annotated[ResultsFileManagement,
        Depends(result_file_management_dependency)]) -> PredictUmolDockingFeature:
    return PredictUmolDockingFeature(target_file_management=target_file_management,
                                     ligand_file_management=ligand_file_management,
                                     result_file_management=result_file_management,
                                     settings=settings)


def update_docking_params_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> UpdateDockingParamsFeature:
    return UpdateDockingParamsFeature(file_management=result_file_management)


def get_diffdock_params_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> GetDiffDockParamsFeature:
    return GetDiffDockParamsFeature(file_management=result_file_management)


def update_diffdock_params_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> UpdateDiffDockParamsFeature:
    return UpdateDiffDockParamsFeature(file_management=result_file_management)


def predict_diffdock_docking_dependency(
        settings: Annotated[Settings, Depends(settings_dependency)],
        target_file_management: Annotated[TargetsFileManagement,
        Depends(target_file_management_dependency)], ligand_file_management: Annotated[LigandsFileManagement,
        Depends(ligand_file_management_dependency)], result_file_management: Annotated[ResultsFileManagement,
        Depends(result_file_management_dependency)]) -> PredictDiffDockDockingFeature:
    return PredictDiffDockDockingFeature(target_file_management=target_file_management,
                                         ligand_file_management=ligand_file_management,
                                         result_file_management=result_file_management,
                                         settings=settings)


def get_umol_docking_result_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> GetUmolDockingResultsFeature:
    return GetUmolDockingResultsFeature(result_file_management=result_file_management)


def get_diffdock_docking_result_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> GetDiffDockDockingResultsFeature:
    return GetDiffDockDockingResultsFeature(result_file_management=result_file_management)


def get_diffdock_ligand_sdf_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]) -> GetDiffDockLigandSdfFeature:
    return GetDiffDockLigandSdfFeature(result_file_management=result_file_management)


def create_experiment_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> CreateExperimentFeature:
    return CreateExperimentFeature(file_management=file_management)
