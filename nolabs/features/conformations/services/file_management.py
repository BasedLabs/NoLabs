import dataclasses
import datetime
import glob
import json
import os.path
import shutil
from typing import Dict

import nolabs
from nolabs.api_models.conformations import RunSimulationsRequest

from nolabs.domain.experiment import ExperimentId, ExperimentName, ExperimentMetadata
from nolabs.infrastructure.settings import Settings
from nolabs.utils.datetime_utils import DateTimeUtils


class FileManagement:
    def __init__(self, settings: Settings, dt_utils: DateTimeUtils):
        self._settings = settings
        self._dt_utils = dt_utils
        self._simulations_file_name = 'simulations.pdb'

    def ensure_experiments_folder_exists(self):
        if not os.path.isdir(self._settings.conformations_experiments_folder):
            os.mkdir(self._settings.conformations_experiments_folder)

    def experiment_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.conformations_experiments_folder, experiment_id.value)

    def create_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        if not os.path.isdir(experiment_folder):
            os.mkdir(experiment_folder)

    def delete_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        shutil.rmtree(experiment_folder, ignore_errors=True)

    def get_experiment_metadata(self) -> ExperimentMetadata:
        metadata_file = os.path.join(self._settings.conformations_experiments_folder,
                                     self._settings.conformations_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        id = ExperimentId(metadata['id'])
        return ExperimentMetadata(
            id=id,
            name=ExperimentName(metadata['name']),
            date=metadata['date'],
            properties=metadata['properties']
        )

    def get_all_experiments_metadata(self) -> Dict[ExperimentId, ExperimentMetadata]:
        metadata_files = os.path.join(self._settings.conformations_experiments_folder,
                                      f'*{self._settings.conformations_metadata_file_name}')
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

    def update_metadata(self, experiment_id: ExperimentId, experiment_name: ExperimentName,
                        run_simulations_request: RunSimulationsRequest):
        j = {
            'id': experiment_id,
            'name': experiment_name,
            'date': self._dt_utils.utcnow(),
            'properties': dataclasses.asdict(run_simulations_request)
        }

        # Remove pdb file from the metadata
        del j['parameters']['pdbFile']  # type: ignore

        metadata_file_path = os.path.join(self.experiment_folder(experiment_id),
                                          self._settings.conformations_metadata_file_name)
        with open(metadata_file_path, 'w', encoding='utf-8') as f:
            json.dump(j, f, ensure_ascii=False, indent=4)

    def save_experiment(self, experiment_id: ExperimentId, pdb_content: str):
        experiment_folder = self.experiment_folder(experiment_id)
        results_pdb_path = os.path.join(experiment_folder, self._simulations_file_name)
        with open(results_pdb_path, 'w', encoding='utf-8') as simulations_f:
            simulations_f.write(pdb_content)

    def get_experiment_data(self, experiment_id: ExperimentId) -> str:
        experiment_folder = self.experiment_folder(experiment_id)
        result_pdb_path = os.path.join(experiment_folder, self._simulations_file_name)
        with open(result_pdb_path, 'r', encoding='utf-8') as simulations_f:
            return simulations_f.read()

    def change_experiment_name(self, experiment_id: ExperimentId, experiment_name: ExperimentName):
        metadata_file = os.path.join(self._settings.conformations_experiments_folder,
                                     self._settings.conformations_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        metadata['name'] = experiment_name.value
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, ensure_ascii=False, indent=4)

