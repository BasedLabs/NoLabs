from nolabs.api_models.drug_discovery import UploadTargetLigandRequest, UploadTargetLigandResponse, \
    GetTargetLigandsListRequest, GetTargetLigandsListResponse, DeleteTargetLigandRequest, DeleteTargetLigandResponse, \
    GetTargetLigandDataRequest, GetTargetLigandDataResponse, GetTargetLigandMetaDataRequest, GetTargetLigandMetaDataResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement


class UploadTargetLigandFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: UploadTargetLigandRequest) -> UploadTargetLigandResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_file = request.sdf_file

        ligand_metadata = self._file_management.store_target_ligand(experiment_id, target_id, ligand_file)

        return UploadTargetLigandResponse(ligand_meta_data=ligand_metadata)


class DeleteTargetLigandFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: DeleteTargetLigandRequest) -> DeleteTargetLigandResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        deleted_ligand_id = self._file_management.delete_target_ligand(experiment_id, target_id, ligand_id)

        return DeleteTargetLigandResponse(ligand_id=deleted_ligand_id.value)


class GetTargetLigandsListFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetLigandsListRequest) -> GetTargetLigandsListResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        ligands = self._file_management.get_target_ligands_list(experiment_id, target_id)

        return GetTargetLigandsListResponse(ligands=ligands)


class GetTargetLigandMetaDataFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetLigandMetaDataRequest) -> GetTargetLigandMetaDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        ligand_metadata = self._file_management.get_target_ligand_metadata(experiment_id, target_id, ligand_id)

        return GetTargetLigandMetaDataResponse(ligand_id=ligand_id.value,
                                         ligand_name=ligand_metadata.ligand_name)


class GetTargetLigandDataFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetTargetLigandDataRequest) -> GetTargetLigandDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        ligand_name, sdf_contents, ligand_smiles = self._file_management.get_target_ligand_data(experiment_id, target_id,
                                                                                         ligand_id)

        return GetTargetLigandDataResponse(ligand_id=ligand_id.value,
                                     ligand_name=ligand_name,
                                     ligand_sdf=sdf_contents,
                                     ligand_smiles=ligand_smiles)
