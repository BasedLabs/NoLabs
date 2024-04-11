from nolabs.api_models.experiment import ChangeExperimentNameRequest
from nolabs.modules.file_management_base import ExperimentsFileManagementBase
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.exceptions import NoLabsException, ErrorCodes


class ChangeExperimentNameFeature:
    def __init__(self, file_management: ExperimentsFileManagementBase):
        self._file_management = file_management

    def handle(self, request: ChangeExperimentNameRequest):
        assert request

        experiment_id = ExperimentId(request.id)
        experiment_name = ExperimentName(request.name)

        if not self._file_management.metadata_exists(experiment_id):
            raise NoLabsException(messages=["Experiment does not exist"], error_code=ErrorCodes.experiment_not_found)

        self._file_management.change_experiment_name(experiment_id, experiment_name)