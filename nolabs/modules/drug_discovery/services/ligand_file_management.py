import json
import os.path
import shutil
from typing import List, Dict

from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.services.file_management import FileManagement
from nolabs.infrastructure.settings import Settings
from nolabs.utils.generate_2d_drug import generate_png_from_smiles, image_file_to_base64

from nolabs.utils.sdf import SDFReader, SDFWriter
from nolabs.utils.uuid_utils import generate_uuid
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.api_models.drug_discovery import LigandMetaData
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement

from fastapi import UploadFile


class LigandsFileManagement:
    """
    Each ligand is attached to a certain Target, so LigandsFileManagement is dependent on TargetsFilemanagement
    """
    def __init__(self, settings: Settings, experiments_file_management: FileManagement, targets_file_management: TargetsFileManagement):
        self._settings = settings
        self.sdf_reader = SDFReader()
        self.sdf_writer = SDFWriter()
        self.experiments_file_management = experiments_file_management
        self._targets_file_management = targets_file_management

    def ensure_target_ligands_folder_exists(self, experiment_id: ExperimentId, target_id: TargetId) -> None:
        self.experiments_file_management.ensure_experiment_folder_exists(experiment_id)
        if not os.path.isdir(self.target_ligands_folder(experiment_id, target_id)):
            os.mkdir(self.target_ligands_folder(experiment_id, target_id))

    def ensure_lone_ligands_folder_exists(self, experiment_id: ExperimentId) -> None:
        if not os.path.exists(self.lone_ligands_folder(experiment_id)):
            os.mkdir(self.lone_ligands_folder(experiment_id))

    def ensure_target_ligand_folder_exists(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId):
        ligand_folder = self.target_ligand_folder(experiment_id, target_id, ligand_id)
        if not os.path.isdir(ligand_folder):
            os.mkdir(ligand_folder)

    def ensure_lone_ligand_folder_exists(self, experiment_id: ExperimentId, ligand_id: LigandId):
        ligand_folder = self.lone_ligand_folder(experiment_id, ligand_id)
        if not os.path.isdir(ligand_folder):
            os.mkdir(ligand_folder)

    def lone_ligands_folder(self, experiment_id: ExperimentId) -> str:
        experiments_folder = self.experiments_file_management.experiment_folder(experiment_id)
        return os.path.join(experiments_folder, 'ligands')

    def lone_ligand_folder(self, experiment_id: ExperimentId, ligand_id: LigandId) -> str:
        return os.path.join(self.lone_ligands_folder(experiment_id), ligand_id.value)

    def target_ligands_folder(self, experiment_id: ExperimentId, target_id: TargetId) -> str:
        target_folder = self._targets_file_management.target_folder(experiment_id, target_id)
        return os.path.join(target_folder, 'ligands')

    def target_ligand_folder(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId) -> str:
        return os.path.join(self.target_ligands_folder(experiment_id, target_id), ligand_id.value)

    def create_target_ligand_folder(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId):
        self.ensure_target_ligands_folder_exists(experiment_id, target_id)
        ligands_dir = self.target_ligands_folder(experiment_id, target_id)
        os.mkdir(os.path.join(ligands_dir, ligand_id.value))

    def create_lone_ligand_folder(self, experiment_id: ExperimentId, ligand_id: LigandId):
        self.ensure_lone_ligands_folder_exists(experiment_id)
        ligands_dir = self.lone_ligands_folder(experiment_id)
        os.mkdir(os.path.join(ligands_dir, ligand_id.value))

    def store_target_ligand(self, experiment_id: ExperimentId, target_id: TargetId,
                            sdf_file: UploadFile) -> LigandMetaData:
        self.ensure_target_ligands_folder_exists(experiment_id, target_id)

        original_filename = sdf_file.filename
        contents = sdf_file.file.read().decode('utf-8')

        ligand_id = LigandId(generate_uuid())
        ligand_name = os.path.splitext(os.path.basename(original_filename))[0]

        self.ensure_target_ligand_folder_exists(experiment_id, target_id, ligand_id)
        ligand_file = os.path.join(self.target_ligand_folder(experiment_id, target_id, ligand_id), original_filename)
        self.sdf_writer.write_sdf(contents, ligand_file)

        smiles = self.sdf_reader.get_smiles_from_sdf(ligand_file)

        image2D = generate_png_from_smiles(smiles)

        image2D.save(os.path.join(self.target_ligand_folder(experiment_id, target_id, ligand_id), "ligand.png"))

        self.update_target_ligand_metadata(experiment_id, target_id, ligand_id, "id", ligand_id.value)
        self.update_target_ligand_metadata(experiment_id, target_id, ligand_id, "name", ligand_name)

        image = image_file_to_base64(os.path.join(self.target_ligand_folder(experiment_id, target_id, ligand_id), "ligand.png"))

        return LigandMetaData(ligand_id=ligand_id.value, ligand_name=ligand_name, image=image)


    def store_lone_ligand(self, experiment_id: ExperimentId, sdf_file: UploadFile,
                          additional_metadata: Dict[str, str] = None) -> LigandMetaData:
        self.ensure_lone_ligands_folder_exists(experiment_id)

        original_filename = sdf_file.filename
        contents = sdf_file.file.read().decode('utf-8')

        ligand_id = LigandId(generate_uuid())
        ligand_name = os.path.splitext(os.path.basename(original_filename))[0]

        self.ensure_lone_ligand_folder_exists(experiment_id, ligand_id)
        ligand_file = os.path.join(self.lone_ligand_folder(experiment_id, ligand_id), original_filename)
        self.sdf_writer.write_sdf(contents, ligand_file)

        smiles = self.sdf_reader.get_smiles_from_sdf(ligand_file)

        image2D = generate_png_from_smiles(smiles)

        image2D.save(os.path.join(self.lone_ligand_folder(experiment_id, ligand_id), "ligand.png"))

        self.update_lone_ligand_metadata(experiment_id, ligand_id, "id", ligand_id.value)
        self.update_lone_ligand_metadata(experiment_id, ligand_id, "name", ligand_name)

        if additional_metadata:
            for key, value in additional_metadata.items():
                self.update_lone_ligand_metadata(experiment_id, ligand_id, key, value)

        image = image_file_to_base64(os.path.join(self.lone_ligand_folder(experiment_id, ligand_id), "ligand.png"))

        if not additional_metadata or "link" not in additional_metadata:
            return LigandMetaData(ligand_id=ligand_id.value, ligand_name=ligand_name, image=image)
        else:
            return LigandMetaData(ligand_id=ligand_id.value, ligand_name=ligand_name,
                                  link=additional_metadata["link"], image=image)


    def delete_target_ligand(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId):
        self.ensure_target_ligands_folder_exists(experiment_id, target_id)
        self.ensure_target_ligand_folder_exists(experiment_id, target_id, ligand_id)
        ligand_folder = self.target_ligand_folder(experiment_id, target_id, ligand_id)
        shutil.rmtree(ligand_folder)

        return ligand_id

    def delete_lone_ligand(self, experiment_id: ExperimentId, ligand_id: LigandId):
        self.ensure_lone_ligands_folder_exists(experiment_id)
        self.ensure_lone_ligand_folder_exists(experiment_id, ligand_id)
        ligand_folder = self.lone_ligand_folder(experiment_id, ligand_id)
        shutil.rmtree(ligand_folder)

        return ligand_id

    def update_target_ligand_metadata(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId, key: str, value: str):
        metadata_file = os.path.join(self.target_ligand_folder(experiment_id, target_id, ligand_id),
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

    def update_lone_ligand_metadata(self, experiment_id: ExperimentId, ligand_id: LigandId, key: str, value: str):
        metadata_file = os.path.join(self.lone_ligand_folder(experiment_id, ligand_id),
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

    def get_target_ligand_metadata(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId) -> LigandMetaData:
        metadata_file = os.path.join(self.target_ligand_folder(experiment_id, target_id, ligand_id),
                                     self._settings.drug_discovery_ligand_metadata_file_name)
        with open(metadata_file, "r") as f:
            metadata = json.load(f)
        ligand_metadata = LigandMetaData(ligand_id=metadata["id"], ligand_name=metadata["name"])
        return ligand_metadata

    def get_lone_ligand_metadata(self, experiment_id: ExperimentId, ligand_id: LigandId) -> LigandMetaData:
        metadata_file = os.path.join(self.lone_ligand_folder(experiment_id, ligand_id),
                                     self._settings.drug_discovery_ligand_metadata_file_name)
        with open(metadata_file, "r") as f:
            metadata = json.load(f)
        image = image_file_to_base64(os.path.join(self.lone_ligand_folder(experiment_id, ligand_id), "ligand.png"))
        if "link" in metadata:
            ligand_metadata = LigandMetaData(ligand_id=metadata["id"], ligand_name=metadata["name"],
                                             link=metadata["link"], image=image)
            return ligand_metadata
        else:
            ligand_metadata = LigandMetaData(ligand_id=metadata["id"], ligand_name=metadata["name"], image=image)
            return ligand_metadata

    def get_target_ligands_list(self, experiment_id: ExperimentId, target_id: TargetId) -> List[LigandMetaData]:
        self.ensure_target_ligands_folder_exists(experiment_id, target_id)
        ligands_folder = self.target_ligands_folder(experiment_id, target_id)
        result_list = []
        for t_id in os.listdir(ligands_folder):
            if os.path.isdir(os.path.join(ligands_folder, t_id)):
                ligand_id = LigandId(t_id)
                ligand_metadata = self.get_target_ligand_metadata(experiment_id, target_id, ligand_id)
                result_list.append(ligand_metadata)

        return result_list

    def get_lone_ligands_list(self, experiment_id: ExperimentId) -> List[LigandMetaData]:
        self.ensure_lone_ligands_folder_exists(experiment_id)
        ligands_folder = self.lone_ligands_folder(experiment_id)
        result_list = []
        for t_id in os.listdir(ligands_folder):
            if os.path.isdir(os.path.join(ligands_folder, t_id)):
                ligand_id = LigandId(t_id)
                ligand_metadata = self.get_lone_ligand_metadata(experiment_id, ligand_id)
                result_list.append(ligand_metadata)

        return result_list

    def get_target_ligand_data(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId) -> tuple[str, str, str]:
        ligand_folder = self.target_ligand_folder(experiment_id, target_id, ligand_id)
        ligand_metadata = self.get_target_ligand_metadata(experiment_id, target_id, ligand_id)
        ligand_name = ligand_metadata.ligand_name

        sdf_contents = self.sdf_reader.read_sdf(os.path.join(ligand_folder, ligand_name + ".sdf"))
        ligand_smiles = self.sdf_reader.get_smiles_from_sdf(os.path.join(ligand_folder, ligand_name + ".sdf"))

        return ligand_name, sdf_contents, ligand_smiles


    def get_lone_ligand_data(self, experiment_id: ExperimentId, ligand_id: LigandId) -> tuple[str, str, str]:
        ligand_folder = self.lone_ligand_folder(experiment_id, ligand_id)
        ligand_metadata = self.get_lone_ligand_metadata(experiment_id, ligand_id)
        ligand_name = ligand_metadata.ligand_name

        sdf_contents = self.sdf_reader.read_sdf(os.path.join(ligand_folder, ligand_name + ".sdf"))
        ligand_smiles = self.sdf_reader.get_smiles_from_sdf(os.path.join(ligand_folder, ligand_name + ".sdf"))

        return ligand_name, sdf_contents, ligand_smiles