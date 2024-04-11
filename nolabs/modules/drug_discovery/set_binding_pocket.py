from nolabs.api_models.drug_discovery import SetTargetBindingPocketRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement


class SetBindingPocketFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: SetTargetBindingPocketRequest):
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        pocket_ids = request.pocket_ids

        self._file_management.store_binding_pocket(experiment_id, target_id, pocket_ids)