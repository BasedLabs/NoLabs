from nolabs.api_models.drug_discovery import UploadLigandRequest, UploadLigandResponse, \
    GetLigandsListRequest, GetLigandsListResponse, DeleteLigandRequest, DeleteLigandResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.ligand import LigandId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.ligand_file_management import LigandsFileManagement


class UploadLigandFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: UploadLigandRequest) -> UploadLigandResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_file = request.sdf_file

        ligand_metadata = self._file_management.store_ligand(experiment_id, target_id, ligand_file)

        return UploadLigandResponse(ligand_meta_data=ligand_metadata)


class DeleteLigandFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: DeleteLigandRequest) -> DeleteLigandResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        deleted_ligand_id = self._file_management.delete_ligand(experiment_id, target_id, ligand_id)

        return DeleteLigandResponse(ligand_id=deleted_ligand_id.value)


class GetLigandsListFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetLigandsListRequest) -> GetLigandsListResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        ligands = self._file_management.get_ligands_list(experiment_id, target_id)

        return GetLigandsListResponse(ligands=ligands)