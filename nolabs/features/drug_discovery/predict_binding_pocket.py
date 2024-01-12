import p2rank_microservice as microservice

from nolabs.api_models.drug_discovery import PredictBindingPocketRequest, PredictBindingPocketResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.targets_file_management import TargetsFileManagement

class GetBindingPocketFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: PredictBindingPocketRequest) -> PredictBindingPocketResponse:
        assert request

        experiment_id = ExperimentId(request.experimentId)
        target_id = TargetId(request.targetId)

        client = microservice.DefaultApi(api_client=microservice.ApiClient())
        pdb_contents = self._file_management.get_pdb_contents(experiment_id, target_id)
        pdb_fixer_request = microservice.RunP2RankPredictionRequest(pdb_contents=pdb_contents)
        pocket_ids = client.predict(request=pdb_fixer_request).pocketIds

        pocket_ids = self._file_management.store_binding_pocket(experiment_id, target_id, pocket_ids)

        return PredictBindingPocketResponse(pocketIds=pocket_ids)