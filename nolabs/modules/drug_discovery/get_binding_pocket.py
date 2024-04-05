from nolabs.api_models.drug_discovery import GetTargetBindingPocketRequest, GetTargetBindingPocketResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement


class GetBindingPocketFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetBindingPocketRequest) -> GetTargetBindingPocketResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        pocket_ids = self._file_management.get_binding_pocket(experiment_id, target_id)

        return GetTargetBindingPocketResponse(pocket_ids=pocket_ids)
