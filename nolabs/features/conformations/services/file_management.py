import dataclasses
import json
import os.path

from nolabs.api_models.conformations import RunSimulationsRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings
from nolabs.utils import utcnow


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.conformations_experiments_folder, settings.conformations_metadata_file_name)
        self._settings = settings
        self._simulations_file_name = settings.conformations_simulations_file_name
        self.ensure_experiments_folder_exists()

    def update_metadata(self, experiment_id: ExperimentId, experiment_name: ExperimentName,
                        run_simulations_request: RunSimulationsRequest):
        j = {
            'id': experiment_id.value,
            'name': experiment_name.value,
            'date': str(utcnow()),
            'properties': dataclasses.asdict(run_simulations_request)
        }

        # Remove pdb file from the metadata
        del j['properties']['pdb_file']  # type: ignore
        del j['properties']['integrator']  # type: ignore

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





