from nolabs.infrastructure.settings import Settings
from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.conformations import RunSimulationsResponse


class GetResultsFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self, experiment_id: str) -> RunSimulationsResponse:
        assert experiment_id

        experiment_id = ExperimentId(experiment_id)
        simulations_data = self._file_management.get_experiment_data(experiment_id)
        return RunSimulationsResponse(pdbContent=simulations_data, errors=[])

