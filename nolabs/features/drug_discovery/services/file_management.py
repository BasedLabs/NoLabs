import glob
import json
import os.path
import shutil
from typing import Dict

from nolabs.domain.experiment import ExperimentId, ExperimentName, ExperimentMetadata
from nolabs.infrastructure.settings import Settings
from nolabs.utils.datetime_utils import DateTimeUtils


class FileManagement:
    def __init__(self, settings: Settings, dt_utils: DateTimeUtils):
        self._settings = settings
        self._dt_utils = dt_utils

    def ensure_experiments_folder_exists(self):
        if not os.path.isdir(self._settings.drug_discovery_experiments_folder):
            os.mkdir(self._settings.drug_discovery_experiments_folder)

    def ensure_targets_folder_exists(self, experiment_id: ExperimentId):
        if not os.path.isdir(self.targets_folder(experiment_id)):
            os.mkdir(self.targets_folder(experiment_id))

    def ensure_ligands_folder_exists(self, experiment_id: ExperimentId):
        if not os.path.isdir(self.ligands_folder(experiment_id)):
            os.mkdir(self.ligands_folder(experiment_id))

    def ensure_results_folder_exists(self, experiment_id: ExperimentId):
        if not os.path.isdir(self.results_folder(experiment_id)):
            os.mkdir(self.results_folder(experiment_id))

    def experiment_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.drug_discovery_experiments_folder, experiment_id.value)

    def targets_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.drug_discovery_experiments_folder, experiment_id.value, 'targets')

    def ligands_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.drug_discovery_experiments_folder, experiment_id.value, 'ligands')

    def results_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.drug_discovery_experiments_folder, experiment_id.value, 'results')

    def create_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        if not os.path.isdir(experiment_folder):
            os.mkdir(experiment_folder)

    def create_target_folder(self, experiment_id: ExperimentId, target_id: str):
        self.ensure_targets_folder_exists(experiment_id)
        targets_dir = self.targets_folder(experiment_id)
        os.mkdir(os.path.join(targets_dir, target_id))


    def delete_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        shutil.rmtree(experiment_folder, ignore_errors=True)

    def get_experiment_metadata(self) -> ExperimentMetadata:
        metadata_file = os.path.join(self._settings.drug_discovery_experiments_folder,
                                     self._settings.drug_discovery_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        id = ExperimentId(metadata['id'])
        return ExperimentMetadata(
            id=id,
            name=ExperimentName(metadata['name']),
            date=metadata['date'],
            properties=metadata['properties']
        )

    def get_all_experiments_metadata(self) -> Dict[ExperimentId, ExperimentMetadata]:
        metadata_files = os.path.join(self._settings.drug_discovery_experiments_folder,
                                      f'*{self._settings.drug_discovery_metadata_file_name}')
        d = {}
        for metadata_file in glob.glob(metadata_files):
            metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
            id = ExperimentId(metadata['id'])
            d[id] = ExperimentMetadata(
                id=id,
                name=ExperimentName(metadata['name']),
                date=metadata['date'],
                properties=metadata['properties']
            )
        return d

    def update_metadata(self, experiment_id: ExperimentId, experiment_name: ExperimentName):
        j = {
            'id': experiment_id,
            'name': experiment_name,
            'date': self._dt_utils.utcnow(),
            'properties': {}
        }

        metadata_file_path = os.path.join(self.experiment_folder(experiment_id),
                                          self._settings.drug_discovery_metadata_file_name)
        with open(metadata_file_path, 'w', encoding='utf-8') as f:

            json.dump(j, f, ensure_ascii=False, indent=4)

    def change_experiment_name(self, experiment_id: ExperimentId, experiment_name: ExperimentName):
        metadata_file = os.path.join(self._settings.drug_discovery_experiments_folder,
                                     experiment_id.value,
                                     self._settings.drug_discovery_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        metadata['name'] = experiment_name.value
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, ensure_ascii=False, indent=4)

