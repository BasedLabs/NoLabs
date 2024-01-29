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
        super().__init__(settings.drug_discovery_experiments_folder,
                         settings.drug_discovery_experiment_metadata_file_name)

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
        return os.path.join(self._experiments_folder, experiment_id.value)

    def targets_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._experiments_folder, experiment_id.value, 'targets')

    def ligands_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._experiments_folder, experiment_id.value, 'ligands')

    def results_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._experiments_folder, experiment_id.value, 'results')

    def create_target_folder(self, experiment_id: ExperimentId, target_id: str):
        self.ensure_targets_folder_exists(experiment_id)
        targets_dir = self.targets_folder(experiment_id)
        os.mkdir(os.path.join(targets_dir, target_id))
