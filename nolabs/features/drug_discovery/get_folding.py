from nolabs.api_models.drug_discovery import GetFoldingRequest, GetFoldingResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement


class GetFoldedStructureFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetFoldingRequest) -> GetFoldingResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        pdb_contents = self._file_management.get_pdb_contents(experiment_id, target_id)

        return GetFoldingResponse(pdb_contents=pdb_contents)
