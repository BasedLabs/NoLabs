import esmfold_light_microservice as microservice

from esmfold_light_microservice import ApiClient, DefaultApi, Configuration
from mypy.errorcodes import ErrorCode

from nolabs.api_models.drug_discovery import PredictFoldingRequest, PredictFoldingResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings


class PredictEsmFoldLightFeature:
    def __init__(self, file_management: TargetsFileManagement,
                 settings: Settings):
        self._file_management = file_management
        self._settings = settings

    def handle(self, request: PredictFoldingRequest) -> PredictFoldingResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        configuration = Configuration(
            host=self._settings.esmfold_light_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            _, sequence, _ = self._file_management.get_target_data(experiment_id, target_id)
            if len(sequence) > 400:
                raise NoLabsException(messages=["Light folding does not support sequences longer than 400. Please use "
                                                "other folding backend"],
                                      error_code=ErrorCodes.drug_discovery_folding_error)
            else:
                request = microservice.RunEsmFoldPredictionRequest(protein_sequence=sequence)
                pdb_content = api_instance.predict_through_api_run_folding_post(run_esm_fold_prediction_request=request).pdb_content

                self._file_management.store_pdb_contents(experiment_id, target_id, pdb_content)

        return PredictFoldingResponse(pdb_content=pdb_content)
