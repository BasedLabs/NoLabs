from typing import Annotated

from fastapi import Depends

from nolabs.controllers.common_dependencies import settings_dependency
from nolabs.modules.experiment.change_experiment_name import ChangeExperimentNameFeature
from nolabs.modules.experiment.create_experiment import CreateExperimentFeature
from nolabs.modules.experiment.get_experiments import GetExperimentsFeature
from nolabs.modules.small_molecules_design.delete_experiment import DeleteExperimentFeature
from nolabs.modules.small_molecules_design.get_experiment import GetExperimentFeature
from nolabs.modules.small_molecules_design.get_logs import GetLogsFeature
from nolabs.modules.small_molecules_design.services.file_management import FileManagement
from nolabs.modules.small_molecules_design.get_smiles import GetSmilesFeature
from nolabs.modules.small_molecules_design.save_properties_feature import SavePropertiesFeature
from nolabs.modules.small_molecules_design.start_learning import StartLearningExperimentFeature
from nolabs.modules.small_molecules_design.start_sampling import StartSamplingExperimentFeature
from nolabs.modules.small_molecules_design.stop_experiment import StopExperimentFeature
from nolabs.modules.small_molecules_design.experiment_status import GetExperimentStatusFeature
from nolabs.infrastructure.settings import Settings


def file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> FileManagement:
    return FileManagement(settings=settings)


def experiment_status_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                                 file_management: Annotated[
                                     FileManagement, Depends(
                                         file_management_dependency)]) -> GetExperimentStatusFeature:
    return GetExperimentStatusFeature(settings=settings, file_management=file_management)


def start_learning_experiment_dependency(
        settings: Annotated[Settings, Depends(settings_dependency)]) -> StartLearningExperimentFeature:
    return StartLearningExperimentFeature(settings=settings)


def start_sampling_experiment_dependency(
        settings: Annotated[Settings, Depends(settings_dependency)]) -> StartSamplingExperimentFeature:
    return StartSamplingExperimentFeature(settings=settings)


def stop_experiment_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> StopExperimentFeature:
    return StopExperimentFeature(settings=settings)


def save_properties_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                               file_management: Annotated[
                                   FileManagement, Depends(file_management_dependency)]) -> SavePropertiesFeature:
    return SavePropertiesFeature(settings=settings, file_management=file_management)


def get_smiles_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                          file_management: Annotated[
                              FileManagement, Depends(file_management_dependency)]
                          ) -> GetSmilesFeature:
    return GetSmilesFeature(file_management=file_management, settings=settings)


def get_logs_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                        file_management: Annotated[
                            FileManagement, Depends(file_management_dependency)]) -> GetLogsFeature:
    return GetLogsFeature(settings=settings, file_management=file_management)


def get_experiments_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> GetExperimentsFeature:
    return GetExperimentsFeature(file_management=file_management)


def get_experiment_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                              file_management: Annotated[
                                  FileManagement, Depends(file_management_dependency)]) -> GetExperimentFeature:
    return GetExperimentFeature(settings=settings, file_management=file_management)


def delete_experiment_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                                 file_management: Annotated[
                                     FileManagement, Depends(file_management_dependency)]) -> DeleteExperimentFeature:
    return DeleteExperimentFeature(settings=settings, file_management=file_management)


def change_experiment_name_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]
) -> ChangeExperimentNameFeature:
    return ChangeExperimentNameFeature(
        file_management=file_management
    )


def create_experiment_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> CreateExperimentFeature:
    return CreateExperimentFeature(file_management=file_management)
