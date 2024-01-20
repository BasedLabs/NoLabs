from nolabs.api_models.protein_design import ChangeExperimentNameRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.features.protein_design.services.file_management import FileManagement
from nolabs.exceptions import NoLabsException, ErrorCodes


class ChangeExperimentNameFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, request: ChangeExperimentNameRequest):
        assert request

        experiment_id = ExperimentId(request.id)
        experiment_name = ExperimentName(request.name)

        if not self._file_management.experiment_exists(experiment_id):
            raise NoLabsException(message="Experiment does not exist", error_code=ErrorCodes.experiment_id_not_found)

        self._file_management.change_experiment_name(experiment_id, experiment_name)