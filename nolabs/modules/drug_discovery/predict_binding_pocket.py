import p2rank_microservice as microservice

from p2rank_microservice import ApiClient, DefaultApi, Configuration

from nolabs.api_models.drug_discovery import PredictBindingPocketRequest, PredictBindingPocketResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings


class PredictBindingPocketFeature:
    def __init__(self, file_management: TargetsFileManagement,
                 settings: Settings):
        self._file_management = file_management
        self._settings = settings

    def handle(self, request: PredictBindingPocketRequest) -> PredictBindingPocketResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        folding_method = self._file_management.get_target_metadata(experiment_id, target_id).folding_method

        configuration = Configuration(
            host=self._settings.p2rank_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            pdb_contents = self._file_management.get_pdb_contents(experiment_id, target_id, folding_method)
            p2rank_request = microservice.RunP2RankPredictionRequest(pdb_contents=pdb_contents)
            pocket_ids = api_instance.predict_run_p2rank_post(run_p2_rank_prediction_request=p2rank_request).pocket_ids

            self._file_management.store_binding_pocket(experiment_id, target_id, pocket_ids)

        return PredictBindingPocketResponse(pocket_ids=pocket_ids)
