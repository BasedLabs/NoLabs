from nolabs.api_models.drug_discovery import UploadLigandRequest, UploadLigandResponse, \
    GetLigandsListRequest, GetLigandsListResponse, DeleteLigandRequest, DeleteLigandResponse, \
    GetLigandDataRequest, GetLigandDataResponse, GetLigandMetaDataRequest, GetLigandMetaDataResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.ligand import LigandId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.ligand_file_management import LigandsFileManagement


class UploadLoneLigandFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: UploadLigandRequest) -> UploadLigandResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_file = request.sdf_file

        ligand_metadata = self._file_management.store_lone_ligand(experiment_id, ligand_file)

        return UploadLigandResponse(ligand_meta_data=ligand_metadata)


class DeleteTargetLigandFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: DeleteLigandRequest) -> DeleteLigandResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        deleted_ligand_id = self._file_management.delete_target_ligand(experiment_id, target_id, ligand_id)

        return DeleteLigandResponse(ligand_id=deleted_ligand_id.value)


class GetTargetLigandsListFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetLigandsListRequest) -> GetLigandsListResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        ligands = self._file_management.get_target_ligands_list(experiment_id, target_id)

        return GetLigandsListResponse(ligands=ligands)


class GetTargetLigandMetaDataFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetLigandMetaDataRequest) -> GetLigandMetaDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        ligand_metadata = self._file_management.get_target_ligand_metadata(experiment_id, target_id, ligand_id)

        return GetLigandMetaDataResponse(ligand_id=ligand_id.value,
                                         ligand_name=ligand_metadata.ligand_name)


class GetTargetLigandDataFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetLigandDataRequest) -> GetLigandDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        ligand_name, sdf_contents, ligand_smiles = self._file_management.get_target_ligand_data(experiment_id, target_id,
                                                                                         ligand_id)

        return GetLigandDataResponse(ligand_id=ligand_id.value,
                                     ligand_name=ligand_name,
                                     ligand_sdf=sdf_contents,
                                     ligand_smiles=ligand_smiles)
