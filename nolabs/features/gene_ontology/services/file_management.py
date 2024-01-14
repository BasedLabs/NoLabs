import dataclasses
import glob
import jsonpickle
import json
import os.path
import shutil
from typing import Dict

from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.domain.gene_ontology import OboNode
from nolabs.api_models.conformations import RunSimulationsRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName, ExperimentMetadata
from nolabs.infrastructure.settings import Settings
from nolabs.utils import utcnow


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.go_experiments_folder, settings.go_metadata_file_name)
        self._settings = settings
        self._go_file_name = settings.go_file_name
        self.ensure_experiments_folder_exists()

    def update_metadata(self, experiment_id: ExperimentId, experiment_name: ExperimentName,
                        run_simulations_request: RunSimulationsRequest):
        j = {
            'id': experiment_id,
            'name': experiment_name,
            'date': utcnow(),
            'properties': dataclasses.asdict(run_simulations_request)
        }

        metadata_file_path = os.path.join(self.experiment_folder(experiment_id),
                                          self._settings.go_metadata_file_name)
        with open(metadata_file_path, 'w', encoding='utf-8') as f:

            json.dump(j, f, ensure_ascii=False, indent=4)

    def save_experiment(self, experiment_id: ExperimentId, data: Dict[str, OboNode]):
        experiment_folder = self.experiment_folder(experiment_id)
        results_pdb_path = os.path.join(experiment_folder, self._go_file_name)
        with open(results_pdb_path, 'w', encoding='utf-8') as go_f:
            j = jsonpickle.dumps(data)
            go_f.write(j)

    def get_experiment_data(self, experiment_id: ExperimentId) -> Dict[str, OboNode]:
        experiment_folder = self.experiment_folder(experiment_id)
        result_pdb_path = os.path.join(experiment_folder, self._go_file_name)
        with open(result_pdb_path, 'r', encoding='utf-8') as simulations_f:
            s = simulations_f.read()
            return jsonpickle.decode(s)
