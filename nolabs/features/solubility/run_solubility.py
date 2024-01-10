from solubility_microservice import (RunSolubilityPredictionRequest,
                                     ApiClient, DefaultApi, Configuration)

from nolabs.domain.solubility import SolubilityProbability
from nolabs.infrastructure.settings import Settings
from nolabs.api_models.solubility import RunSolubilityResponse, RunSolubilityRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.features.solubility.services.file_management import FileManagement


class RunSolubilityFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._settings = settings
        self._file_management = file_management

    async def handle(self, request: RunSolubilityRequest) -> RunSolubilityResponse:
        assert request

        configuration = Configuration(
            host=self._settings.solubility_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            result = api_instance.run_solubility_run_solubility_prediction_post(
                run_solubility_prediction_request=RunSolubilityPredictionRequest(
                    protein_sequence=request.aminoAcidSequence
                )
            )
            self._file_management.save_experiment(
                experiment_id=ExperimentId(request.experimentId),
                soluble_probability=SolubilityProbability(result.solubleProbability)
            )
            return RunSolubilityResponse(result.solubleProbability)



