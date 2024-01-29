from typing import Annotated

from fastapi import Depends

from nolabs.controllers.common_dependencies import settings_dependency
from nolabs.features.experiment.create_experiment import CreateExperimentFeature
from nolabs.features.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.features.experiment.change_experiment_name import ChangeExperimentNameFeature
from nolabs.features.drug_discovery.generate_msa import GenerateMsaFeature
from nolabs.features.drug_discovery.services.file_management import FileManagement
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.features.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.features.drug_discovery.services.result_file_management import ResultsFileManagement
from nolabs.features.drug_discovery.target_management import UploadTargetFeature, DeleteTargetFeature, \
    GetTargetsListFeature, GetTargetDataFeature
from nolabs.features.drug_discovery.ligand_management import UploadLigandFeature, DeleteLigandFeature, \
    GetLigandsListFeature, GetLigandDataFeature
from nolabs.features.drug_discovery.get_binding_pocket import GetBindingPocketFeature
from nolabs.features.drug_discovery.predict_binding_pocket import PredictBindingPocketFeature
from nolabs.features.drug_discovery.predict_light_folding import PredictFoldingFeature
from nolabs.features.drug_discovery.get_folding import GetFoldedStructureFeature
from nolabs.features.drug_discovery.register_docking_job import RegisterDockingJobFeature
from nolabs.features.drug_discovery.predict_docking import PredictDockingFeature
from nolabs.features.experiment.get_experiments import GetExperimentsFeature
from nolabs.infrastructure.settings import Settings


def file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> FileManagement:
    return FileManagement(settings=settings)


def add_experiment_feature_dependency(file_management: Annotated[FileManagement, Depends(file_management_dependency)]):
    return AddExperimentFeature(file_management=file_management)


def get_experiments_feature_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> GetExperimentsFeature:
    return GetExperimentsFeature(file_management=file_management)


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
                                      target_file_management: Annotated[TargetsFileManagement,
                                      Depends(target_file_management_dependency)]) -> LigandsFileManagement:
    return LigandsFileManagement(settings=settings, targets_file_management=target_file_management)


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


def get_target_data_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> GetTargetDataFeature:
    return GetTargetDataFeature(file_management=target_file_management)


def get_binding_pocket_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> GetBindingPocketFeature:
    return GetBindingPocketFeature(file_management=target_file_management)


def predict_binding_pocket_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                                      target_file_management: Annotated[TargetsFileManagement,
                                      Depends(target_file_management_dependency)]) -> PredictBindingPocketFeature:
    return PredictBindingPocketFeature(file_management=target_file_management, settings=settings)


def get_folded_structure_dependency(target_file_management: Annotated[TargetsFileManagement,
Depends(target_file_management_dependency)]) -> GetFoldedStructureFeature:
    return GetFoldedStructureFeature(file_management=target_file_management)


def predict_folding_dependency(target_file_management: Annotated[TargetsFileManagement,
                               Depends(target_file_management_dependency)],
                               settings: Annotated[Settings,
                               Depends(settings_dependency)]) -> PredictFoldingFeature:
    return PredictFoldingFeature(file_management=target_file_management, settings=settings)


def upload_ligand_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> UploadLigandFeature:
    return UploadLigandFeature(file_management=ligand_file_management)


def delete_ligand_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> DeleteLigandFeature:
    return DeleteLigandFeature(file_management=ligand_file_management)


def get_ligands_list_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> GetLigandsListFeature:
    return GetLigandsListFeature(file_management=ligand_file_management)


def get_ligand_data_dependency(ligand_file_management: Annotated[LigandsFileManagement,
Depends(ligand_file_management_dependency)]) -> GetLigandDataFeature:
    return GetLigandDataFeature(file_management=ligand_file_management)


def register_docking_job_dependency(result_file_management: Annotated[ResultsFileManagement,
Depends(result_file_management_dependency)]):
    return RegisterDockingJobFeature(file_management=result_file_management)


def predict_docking_dependency(
        settings: Annotated[Settings, Depends(settings_dependency)],
        target_file_management: Annotated[TargetsFileManagement,
        Depends(target_file_management_dependency)], ligand_file_management: Annotated[LigandsFileManagement,
        Depends(ligand_file_management_dependency)], result_file_management: Annotated[ResultsFileManagement,
        Depends(result_file_management_dependency)]) -> PredictDockingFeature:
    return PredictDockingFeature(target_file_management=target_file_management,
                                 ligand_file_management=ligand_file_management,
                                 result_file_management=result_file_management,
                                 settings=settings)

def create_experiment_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> CreateExperimentFeature:
    return CreateExperimentFeature(file_management=file_management)
