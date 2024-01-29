import glob
import json
import os.path
import shutil
from typing import Dict, List

from nolabs.domain.experiment import ExperimentId, ExperimentName, ExperimentMetadata
from nolabs.features.file_management_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings
from nolabs.utils import utcnow, generate_uuid


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__('ASD', 'ASD') # TODO to change
        self._settings = settings
        self.ensure_experiments_folder_exists()

    def ensure_experiments_folder_exists(self):
        if not os.path.isdir(self._settings.drug_discovery_experiments_folder):
            os.makedirs(self._settings.drug_discovery_experiments_folder, exist_ok=True)

    def ensure_targets_folder_exists(self, experiment_id: ExperimentId):
        if not os.path.isdir(self.targets_folder(experiment_id)):
            os.makedirs(self.targets_folder(experiment_id), exist_ok=True)

    def ensure_ligands_folder_exists(self, experiment_id: ExperimentId):
        if not os.path.isdir(self.ligands_folder(experiment_id)):
            os.makedirs(self.ligands_folder(experiment_id), exist_ok=True)

    def ensure_results_folder_exists(self, experiment_id: ExperimentId):
        if not os.path.isdir(self.results_folder(experiment_id)):
            os.makedirs(self.results_folder(experiment_id), exist_ok=True)

    def experiment_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.drug_discovery_experiments_folder, experiment_id.value)

    def targets_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.drug_discovery_experiments_folder, experiment_id.value, 'targets')

    def ligands_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.drug_discovery_experiments_folder, experiment_id.value, 'ligands')

    def results_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.drug_discovery_experiments_folder, experiment_id.value, 'results')

    def create_experiment_folder(self) -> ExperimentMetadata:
        experiment_id = ExperimentId(generate_uuid())
        experiment_folder = self.experiment_folder(experiment_id)
        if not os.path.isdir(experiment_folder):
            os.mkdir(experiment_folder)

        self.set_metadata(experiment_id, ExperimentName("New Experiment"))
        metadata = self.get_metadata(experiment_id)

        return metadata

    def create_target_folder(self, experiment_id: ExperimentId, target_id: str):
        self.ensure_targets_folder_exists(experiment_id)
        targets_dir = self.targets_folder(experiment_id)
        os.mkdir(os.path.join(targets_dir, target_id))


    def delete_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        shutil.rmtree(experiment_folder, ignore_errors=True)

    def get_metadata(self, experiment_id: ExperimentId) -> ExperimentMetadata:
        metadata_file = os.path.join(self.experiment_folder(experiment_id),
                                     self._settings.drug_discovery_experiment_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        id = ExperimentId(metadata['id'])
        return ExperimentMetadata(
            id=id,
            name=ExperimentName(metadata['name']),
            date=metadata['date'],
            properties=metadata['properties']
        )

    def get_all_experiments_metadata(self) -> List[ExperimentMetadata]:
        experiments_list = []
        for exp_id in os.listdir(self._settings.drug_discovery_experiments_folder):
            experiment_id = ExperimentId(exp_id)
            metadata = self.get_metadata(experiment_id)
            experiments_list.append(metadata)

        return experiments_list

    def change_experiment_name(self, experiment_id: ExperimentId, experiment_name: ExperimentName):
        metadata_file = os.path.join(self.experiment_folder(experiment_id),
                                     self._settings.drug_discovery_experiment_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        metadata['name'] = experiment_name.value
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, ensure_ascii=False, indent=4)

