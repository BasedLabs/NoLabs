from nolabs.api_models.drug_discovery import ExperimentMetadataResponse
from nolabs.features.drug_discovery.services.file_management import FileManagement


class AddExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self) -> ExperimentMetadataResponse:
        experiment_metadata = self._file_management.create_experiment_folder()

        return ExperimentMetadataResponse(experiment_id=experiment_metadata.id.value,
                                          experiment_name=experiment_metadata.name.value,
                                          experiment_date=experiment_metadata.date)
