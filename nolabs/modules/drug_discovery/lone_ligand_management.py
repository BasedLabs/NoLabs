from nolabs.api_models.drug_discovery import UploadLoneLigandRequest, UploadLoneLigandResponse, \
    GetLoneLigandsListRequest, GetLoneLigandsListResponse, DeleteLoneLigandRequest, DeleteLoneLigandResponse, \
    GetLoneLigandDataRequest, GetLoneLigandDataResponse, GetLoneLigandMetaDataRequest, GetLoneLigandMetaDataResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement


class UploadLoneLigandFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: UploadLoneLigandRequest) -> UploadLoneLigandResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        ligand_file = request.sdf_file
        metadata = request.metadata

        ligand_metadata = self._file_management.store_lone_ligand(experiment_id, ligand_file, metadata)

        return UploadLoneLigandResponse(ligand_meta_data=ligand_metadata)


class DeleteLoneLigandFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: DeleteLoneLigandRequest) -> DeleteLoneLigandResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        ligand_id = LigandId(request.ligand_id)

        deleted_ligand_id = self._file_management.delete_lone_ligand(experiment_id, ligand_id)

        return DeleteLoneLigandResponse(ligand_id=deleted_ligand_id.value)


class GetLoneLigandsListFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetLoneLigandsListRequest) -> GetLoneLigandsListResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        ligands = self._file_management.get_lone_ligands_list(experiment_id)

        return GetLoneLigandsListResponse(ligands=ligands)


class GetLoneLigandMetaDataFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetLoneLigandMetaDataRequest) -> GetLoneLigandMetaDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        ligand_id = LigandId(request.ligand_id)

        ligand_metadata = self._file_management.get_lone_ligand_metadata(experiment_id, ligand_id)

        return GetLoneLigandMetaDataResponse(ligand_id=ligand_id.value,
                                             ligand_name=ligand_metadata.ligand_name,
                                             image=ligand_metadata.image)


class GetLoneLigandDataFeature:
    def __init__(self, file_management: LigandsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetLoneLigandDataRequest) -> GetLoneLigandDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        ligand_id = LigandId(request.ligand_id)

        ligand_name, sdf_contents, ligand_smiles = self._file_management.get_lone_ligand_data(experiment_id,
                                                                                              ligand_id)

        return GetLoneLigandDataResponse(ligand_id=ligand_id.value,
                                     ligand_name=ligand_name,
                                     ligand_sdf=sdf_contents,
                                     ligand_smiles=ligand_smiles)
