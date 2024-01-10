import dataclasses
import glob
import json
import os.path
import shutil
from typing import Dict

from nolabs.api_models.solubility import RunSolubilityRequest
from nolabs.domain.solubility import SolubilityProbability
from nolabs.domain.experiment import ExperimentId, ExperimentName, ExperimentMetadata
from nolabs.infrastructure.settings import Settings
from nolabs.utils.datetime_utils import DateTimeUtils


class FileManagement:
    def __init__(self, settings: Settings, dt_utils: DateTimeUtils):
        self._settings = settings
        self._dt_utils = dt_utils
        self._solubility_file_name = settings.solubility_file_name

    def ensure_experiments_folder_exists(self):
        if not os.path.isdir(self._settings.solubility_experiments_folder):
            os.mkdir(self._settings.solubility_experiments_folder)

    def experiment_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.solubility_experiments_folder, experiment_id.value)

    def create_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        if not os.path.isdir(experiment_folder):
            os.mkdir(experiment_folder)

    def delete_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        shutil.rmtree(experiment_folder, ignore_errors=True)

    def get_experiment_metadata(self) -> ExperimentMetadata:
        metadata_file = os.path.join(self._settings.solubility_experiments_folder,
                                     self._settings.solubility_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        id = ExperimentId(metadata['id'])
        return ExperimentMetadata(
            id=id,
            name=ExperimentName(metadata['name']),
            date=metadata['date'],
            properties=metadata['properties']
        )

    def get_all_experiments_metadata(self) -> Dict[ExperimentId, ExperimentMetadata]:
        metadata_files = os.path.join(self._settings.solubility_experiments_folder,
                                      f'*{self._settings.solubility_metadata_file_name}')
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
                        run_solubility_request: RunSolubilityRequest):
        j = {
            'id': experiment_id,
            'name': experiment_name,
            'date': self._dt_utils.utcnow(),
            'properties': dataclasses.asdict(run_solubility_request)
        }

        metadata_file_path = os.path.join(self.experiment_folder(experiment_id),
                                          self._settings.go_metadata_file_name)
        with open(metadata_file_path, 'w', encoding='utf-8') as f:
            json.dump(j, f, ensure_ascii=False, indent=4)

    def save_experiment(self, experiment_id: ExperimentId, soluble_probability: SolubilityProbability):
        experiment_folder = self.experiment_folder(experiment_id)
        results_path = os.path.join(experiment_folder, self._solubility_file_name)
        with open(results_path, 'w', encoding='utf-8') as solubility_f:
            solubility_f.write(str(soluble_probability.value))

    def get_experiment_data(self, experiment_id: ExperimentId) -> SolubilityProbability:
        experiment_folder = self.experiment_folder(experiment_id)
        results_path = os.path.join(experiment_folder, self._solubility_file_name)
        with open(results_path, 'r', encoding='utf-8') as solubility_f:
            s = solubility_f.read()
            return SolubilityProbability(
                value=float(s)
            )

    def change_experiment_name(self, experiment_id: ExperimentId, experiment_name: ExperimentName):
        metadata_file = os.path.join(self._settings.solubility_file_name,
                                     experiment_id.value,
                                     self._settings.solubility_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        metadata['name'] = experiment_name.value
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, ensure_ascii=False, indent=4)
