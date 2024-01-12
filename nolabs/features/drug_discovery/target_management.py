from nolabs.api_models.drug_discovery import UploadTargetRequest, UploadTargetResponse, GetTargetsListRequest, GetTargetsListResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.services.targets_file_management import TargetsFileManagement

class UploadTargetFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: UploadTargetRequest) -> UploadTargetResponse:
        assert request

        experiment_id = ExperimentId(request.experimentId)
        protein_file = request.fasta_file

        response = self._file_management.store_target(experiment_id, protein_file)

        return UploadTargetResponse(result=response)

class GetTargetsListFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetsListRequest) -> GetTargetsListResponse:
        assert request

        experiment_id = ExperimentId(request.experimentId)

        targets = self._file_management.get_targets_list(experiment_id)

        return GetTargetsListResponse(targets=targets)