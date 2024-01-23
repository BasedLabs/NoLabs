from typing import Annotated

from fastapi import Depends

from features.experiment.create_experiment import CreateExperimentFeature
from nolabs.controllers.common_dependencies import settings_dependency
from nolabs.features.protein_design import GetExperimentsMetadataFeature, GetExperimentFeature
from nolabs.features.experiment.delete_experiment import DeleteExperimentFeature
from nolabs.features.experiment.change_experiment_name import ChangeExperimentNameFeature
from nolabs.features.protein_design import RunProteinDesignFeature
from nolabs.features.protein_design.services.file_management import FileManagement
from nolabs.infrastructure.settings import Settings


def file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> FileManagement:
    return FileManagement(settings=settings)


def run_protein_design_feature_dependency(file_management: Annotated[FileManagement, Depends(file_management_dependency)],
                                          settings: Annotated[Settings, Depends(settings_dependency)]) -> RunProteinDesignFeature:
    return RunProteinDesignFeature(file_management=file_management,
                                   settings=settings)


def get_experiments_feature_dependency(
        file_management: Annotated[FileManagement, Depends(file_management_dependency)]) -> GetExperimentsMetadataFeature:
    return GetExperimentsMetadataFeature(file_management=file_management)


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

def create_experiment_dependency() -> CreateExperimentFeature:
    return CreateExperimentFeature()
