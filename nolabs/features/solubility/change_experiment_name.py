from nolabs.api_models.solubility import ChangeExperimentNameRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.features.solubility.services.file_management import FileManagement


class ChangeExperimentNameFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, request: ChangeExperimentNameRequest):
        assert request

        experiment_id = ExperimentId(request.id)
        experiment_name = ExperimentName(request.name)

        self._file_management.change_experiment_name(experiment_id, experiment_name)