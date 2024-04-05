from nolabs.modules.file_management_base import ExperimentsFileManagementBase
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes


class DeleteExperimentFeature:
    def __init__(self, file_management: ExperimentsFileManagementBase):
        self._file_management = file_management

    def handle(self, id: str):
        assert id

        experiment_id = ExperimentId(id)

        if not self._file_management.metadata_exists(experiment_id):
            raise NoLabsException(messages=["This experiment not found"], error_code=ErrorCodes.experiment_not_found)

        self._file_management.delete_experiment_folder(experiment_id)