import json
import os.path
import shutil
from typing import List

from nolabs.domain.experiment import ExperimentId
from nolabs.infrastructure.settings import Settings

from nolabs.utils.sdf import SDFReader, SDFWriter
from nolabs.utils.uuid_utils import generate_uuid
from nolabs.features.drug_discovery.data_models.ligand import LigandId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.api_models.drug_discovery import LigandMetaData
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement

from fastapi import UploadFile


class LigandsFileManagement:
    """
    Each ligand is attached to a certain Target, so LigandsFileManagement is dependent on TargetsFilemanagement
    """
    def __init__(self, settings: Settings, targets_file_management: TargetsFileManagement):
        self._settings = settings
        self.sdf_reader = SDFReader()
        self.sdf_writer = SDFWriter()
        self._targets_file_management = targets_file_management

    def ensure_ligands_folder_exists(self, experiment_id: ExperimentId, target_id: TargetId) -> None:
        if not os.path.isdir(self.ligands_folder(experiment_id, target_id)):
            os.mkdir(self.ligands_folder(experiment_id, target_id))

    def ensure_ligand_folder_exists(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId):
        ligand_folder = self.ligand_folder(experiment_id, target_id, ligand_id)
        if not os.path.isdir(ligand_folder):
            os.mkdir(ligand_folder)

    def ligands_folder(self, experiment_id: ExperimentId, target_id: TargetId) -> str:
        target_folder = self._targets_file_management.target_folder(experiment_id, target_id)
        return os.path.join(target_folder, 'ligands')

    def ligand_folder(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId) -> str:
        return os.path.join(self.ligands_folder(experiment_id, target_id), ligand_id.value)

    def create_ligand_folder(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId):
        self.ensure_ligands_folder_exists(experiment_id, target_id)
        ligands_dir = self.ligands_folder(experiment_id, target_id)
        os.mkdir(os.path.join(ligands_dir, ligand_id.value))

    def store_ligand(self, experiment_id: ExperimentId, target_id: TargetId, sdf_file: UploadFile) -> LigandMetaData:
        self.ensure_ligands_folder_exists(experiment_id, target_id)

        original_filename = sdf_file.filename
        contents = sdf_file.file.read().decode('utf-8')

        ligand_id = LigandId(generate_uuid())
        ligand_name = os.path.splitext(os.path.basename(original_filename))[0]

        self.ensure_ligand_folder_exists(experiment_id, target_id, ligand_id)
        ligand_file = os.path.join(self.ligand_folder(experiment_id, target_id, ligand_id), original_filename)
        self.sdf_writer.write_sdf(contents, ligand_file)

        self.update_ligand_metadata(experiment_id, target_id, ligand_id, "id", ligand_id.value)
        self.update_ligand_metadata(experiment_id, target_id, ligand_id, "name", ligand_name)

        return LigandMetaData(ligand_id=ligand_id.value, ligand_name=ligand_name)

    def delete_ligand(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId):
        self.ensure_ligands_folder_exists(experiment_id, target_id)
        self.ensure_ligand_folder_exists(experiment_id, target_id, ligand_id)
        ligand_folder = self.ligand_folder(experiment_id, target_id, ligand_id)
        shutil.rmtree(ligand_folder)

        return ligand_id

    def update_ligand_metadata(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId, key: str, value: str):
        metadata_file = os.path.join(self.ligand_folder(experiment_id, target_id, ligand_id),
                                     self._settings.drug_discovery_ligand_metadata_file_name)
        if os.path.exists(metadata_file):
            with open(metadata_file, "r") as f:
                metadata = json.load(f)
                metadata[key] = value
            with open(metadata_file, "w") as f:
                json.dump(metadata, f)
        else:
            metadata = {key: value}
            with open(metadata_file, "w") as f:
                json.dump(metadata, f)

    def get_ligand_metadata(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId) -> LigandMetaData:
        metadata_file = os.path.join(self.ligand_folder(experiment_id, target_id, ligand_id),
                                     self._settings.drug_discovery_ligand_metadata_file_name)
        with open(metadata_file, "r") as f:
            metadata = json.load(f)
        ligand_metadata = LigandMetaData(ligand_id=metadata["id"], ligand_name=metadata["name"])
        return ligand_metadata

    def get_ligands_list(self, experiment_id: ExperimentId, target_id: TargetId) -> List[LigandMetaData]:
        self.ensure_ligands_folder_exists(experiment_id, target_id)
        ligands_folder = self.ligands_folder(experiment_id, target_id)
        result_list = []
        for t_id in os.listdir(ligands_folder):
            if os.path.isdir(os.path.join(ligands_folder, t_id)):
                ligand_id = LigandId(t_id)
                ligand_metadata = self.get_ligand_metadata(experiment_id, target_id, ligand_id)
                result_list.append(ligand_metadata)

        return result_list

    def get_ligand_data(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId) -> tuple[str, str, str]:
        ligand_folder = self.ligand_folder(experiment_id, target_id, ligand_id)
        ligand_metadata = self.get_ligand_metadata(experiment_id, target_id, ligand_id)
        ligand_name = ligand_metadata.ligand_name

        sdf_contents = self.sdf_reader.read_sdf(os.path.join(ligand_folder, ligand_name + ".sdf"))
        ligand_smiles = self.sdf_reader.get_smiles_from_sdf(os.path.join(ligand_folder, ligand_name + ".sdf"))

        return ligand_name, sdf_contents, ligand_smiles
