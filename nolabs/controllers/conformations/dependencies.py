from typing import Annotated

from fastapi import Depends

from nolabs.features.events_queue import EventsQueue
from nolabs.controllers.common_dependencies import settings_dependency, events_queue_dependency
from nolabs.features.conformations import GetExperimentsFeature, GetExperimentFeature, DeleteExperimentFeature, \
    ChangeExperimentNameFeature
from nolabs.features.conformations import RunSimulationsFeature
from nolabs.features.conformations.services.file_management import FileManagement
from nolabs.infrastructure.settings import Settings


def file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> FileManagement:
    return FileManagement(settings=settings)


def run_simulations_feature_dependency(file_management: Annotated[FileManagement, Depends(file_management_dependency)],
                                       settings: Annotated[Settings, Depends(settings_dependency)],
                                       events_queue: Annotated[EventsQueue, Depends(events_queue_dependency)]
                                       ) -> RunSimulationsFeature:
    return RunSimulationsFeature(file_management=file_management, settings=settings, events_queue=events_queue)


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
