from nolabs.api_models.drug_discovery import UploadTargetRequest, UploadTargetResponse, \
    GetTargetsListRequest, GetTargetsListResponse, DeleteTargetRequest, DeleteTargetResponse, \
    GetTargetDataRequest, GetTargetDataResponse, GetTargetMetaDataRequest, GetTargetMetaDataResponse, \
    UpdateTargetNameRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement


class UploadTargetFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: UploadTargetRequest) -> UploadTargetResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        protein_file = request.fasta_file
        metadata = request.metadata

        response = self._file_management.store_target(experiment_id, protein_file, metadata)

        return UploadTargetResponse(result=response)


class DeleteTargetFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: DeleteTargetRequest) -> DeleteTargetResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        deleted_target_id = self._file_management.delete_target(experiment_id, target_id)

        return DeleteTargetResponse(target_id=deleted_target_id.value)


class GetTargetsListFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetsListRequest) -> GetTargetsListResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        targets = self._file_management.get_targets_list(experiment_id)

        return GetTargetsListResponse(targets=targets)


class GetTargetMetaDataFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetMetaDataRequest) -> GetTargetMetaDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        target_metadata = self._file_management.get_target_metadata(experiment_id, target_id)

        return GetTargetMetaDataResponse(target_id=target_metadata.target_id, target_name=target_metadata.target_name)


class UpdateTargetNameFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: UpdateTargetNameRequest) -> GetTargetMetaDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        target_name = request.target_name

        target_metadata = self._file_management.update_target_name(experiment_id, target_id, target_name)

        return GetTargetMetaDataResponse(target_id=target_metadata.target_id, target_name=target_metadata.target_name)


class GetTargetDataFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetDataRequest) -> GetTargetDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        target_name, sequence, pdb_content = self._file_management.get_target_data(experiment_id, target_id)

        return GetTargetDataResponse(protein_sequence=sequence, protein_pdb=pdb_content)
