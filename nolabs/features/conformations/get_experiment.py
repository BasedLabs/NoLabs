from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.conformations import GetExperimentResponse, ExperimentMetadataResponse
from nolabs.features.conformations.services.file_management import FileManagement
from nolabs.exceptions import NoLabsException, ErrorCodes


class GetExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> GetExperimentResponse:
        assert id

        experiment_id = ExperimentId(id)

        if not self._file_management.experiment_exists(experiment_id):
            raise NoLabsException(message="Experiment does not exist", error_code=ErrorCodes.experiment_id_not_found)

        metadata = self._file_management.get_experiment_metadata(experiment_id)
        data = self._file_management.get_experiment_data(experiment_id)
        return GetExperimentResponse(
            metadata=ExperimentMetadataResponse(
                id=metadata.id.value,
                name=metadata.name.value,
                date=metadata.date
            ),
            data=data
        )
