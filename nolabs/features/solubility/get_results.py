from nolabs.features.solubility.services.file_management import FileManagement
from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.solubility import RunSolubilityResponse


class GetResultsFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> RunSolubilityResponse:
        assert id

        experiment_id = ExperimentId(id)
        solubility_probability = self._file_management.get_experiment_data(experiment_id)

        return RunSolubilityResponse(
            errors=[],
            solubleProbability=solubility_probability.value
        )
