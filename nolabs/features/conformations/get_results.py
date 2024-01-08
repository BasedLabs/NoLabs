from nolabs.features.conformations.services.file_management import FileManagement
from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.conformations import RunSimulationsResponse


class GetResultsFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> RunSimulationsResponse:
        assert id

        experiment_id = ExperimentId(id)
        simulations_data = self._file_management.get_experiment_data(experiment_id)
        return RunSimulationsResponse(pdbContent=simulations_data, errors=[])

