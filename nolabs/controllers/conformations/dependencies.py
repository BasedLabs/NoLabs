from typing import Annotated

from fastapi import Depends

from nolabs.modules.conformations.get_experiment import GetExperimentFeature
from nolabs.modules.conformations.run_simulations import RunSimulationsFeature
from nolabs.modules.experiment.change_experiment_name import ChangeExperimentNameFeature
from nolabs.controllers.common_dependencies import settings_dependency
from nolabs.modules.experiment.create_experiment import CreateExperimentFeature
from nolabs.modules.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.modules.conformations.services.file_management import FileManagement
from nolabs.modules.experiment.get_experiments import GetExperimentsFeature
from nolabs.infrastructure.settings import Settings


def file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> FileManagement:
    return FileManagement(settings=settings)


def run_simulations_feature_dependency(file_management: Annotated[FileManagement, Depends(file_management_dependency)],
                                       settings: Annotated[Settings, Depends(settings_dependency)]
                                       ) -> RunSimulationsFeature:
    return RunSimulationsFeature(file_management=file_management, settings=settings)


def get_experiments_feature_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> GetExperimentsFeature:
    return GetExperimentsFeature(file_management=file_management)


def get_experiment_feature_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> GetExperimentFeature:
    return GetExperimentFeature(file_management=file_management)


def delete_experiment_feature_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> DeleteExperimentFeature:
    return DeleteExperimentFeature(file_management=file_management)


def change_experiment_name_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]
) -> ChangeExperimentNameFeature:
    return ChangeExperimentNameFeature(
        file_management=file_management
    )

def create_experiment_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> CreateExperimentFeature:
    return CreateExperimentFeature(file_management=file_management)
