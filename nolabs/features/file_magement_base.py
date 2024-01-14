import glob
import json
import os
import shutil
from typing import Dict, List

from nolabs.domain.experiment import ExperimentId, ExperimentMetadata, ExperimentName


class ExperimentsFileManagementBase:
    def __init__(self, experiments_folder: str, metadata_file: str):
        self._experiments_folder = experiments_folder
        self._metadata_file = metadata_file

    def ensure_experiments_folder_exists(self):
        if not os.path.isdir(self._experiments_folder):
            os.makedirs(self._experiments_folder, exist_ok=True)

    def ensure_experiment_folder_exists(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        if not os.path.isdir(experiment_folder):
            os.makedirs(experiment_folder, exist_ok=True)

    def experiment_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._experiments_folder, experiment_id.value)

    def create_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        if not os.path.isdir(experiment_folder):
            os.mkdir(experiment_folder)

    def delete_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        shutil.rmtree(experiment_folder, ignore_errors=True)

    def experiment_exists(self, experiment_id: ExperimentId) -> bool:
        return os.path.exists(os.path.join(self.experiment_folder(experiment_id)))

    def get_experiment_metadata(self, experiment_id: ExperimentId) -> ExperimentMetadata:
        metadata_file = os.path.join(self._experiments_folder,
                                     experiment_id.value,
                                     self._metadata_file)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        id = ExperimentId(metadata['id'])
        return ExperimentMetadata(
            id=id,
            name=ExperimentName(metadata['name']),
            date=metadata['date'],
            properties=metadata['properties']
        )

    def change_experiment_name(self, experiment_id: ExperimentId, experiment_name: ExperimentName):
        metadata_file = os.path.join(self.experiment_folder(experiment_id),
                                     self._metadata_file)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        metadata['name'] = experiment_name.value
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, ensure_ascii=False, indent=4)

    def get_all_experiments_metadata(self) -> List[ExperimentMetadata]:
        metadata_files = os.path.join(self._experiments_folder,
                                      '*',
                                      f'*{self._metadata_file}')
        d = []
        for metadata_file in glob.glob(metadata_files):
            metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
            id = ExperimentId(metadata['id'])
            d.append(ExperimentMetadata(
                id=id,
                name=ExperimentName(metadata['name']),
                date=metadata['date'],
                properties=metadata['properties']
            ))
        return d
