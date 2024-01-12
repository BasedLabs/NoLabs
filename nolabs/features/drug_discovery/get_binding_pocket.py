from nolabs.api_models.drug_discovery import GetTargetBindingPocketRequest, GetTargetBindingPocketResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.targets_file_management import TargetsFileManagement

class GetBindingPocketFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetBindingPocketRequest) -> GetTargetBindingPocketResponse:
        assert request

        experiment_id = ExperimentId(request.experimentId)
        target_id = TargetId(request.targetId)

        pocket_ids = self._file_management.get_binding_pocket(experiment_id, target_id)

        return GetTargetBindingPocketResponse(pocketIds=pocket_ids)