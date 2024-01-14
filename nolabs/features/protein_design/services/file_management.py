import dataclasses
import glob
import json
import os.path
import shutil
from typing import Dict, List

from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.api_models.protein_design import RunProteinDesignRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName, ExperimentMetadata
from nolabs.infrastructure.settings import Settings
from nolabs.utils import utcnow


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings, experiments_folder: str, metadata_file: str):
        super().__init__(settings.protein_design_experiments_folder, settings.protein_design_metadata_file_name)
        self._settings = settings
        self._protein_design_file_template = settings.protein_design_file_name
        self._results_folder = 'results'
        self.ensure_experiments_folder_exists()

    def update_metadata(self, experiment_id: ExperimentId, experiment_name: ExperimentName,
                        run_protein_designs_request: RunProteinDesignRequest):
        j = {
            'id': experiment_id,
            'name': experiment_name,
            'date': utcnow(),
            'properties': dataclasses.asdict(run_protein_designs_request)
        }

        # Remove pdb file from the metadata
        del j['parameters']['pdbFile']  # type: ignore

        metadata_file_path = os.path.join(self.experiment_folder(experiment_id),
                                          self._settings.protein_design_metadata_file_name)
        with open(metadata_file_path, 'w', encoding='utf-8') as f:
            json.dump(j, f, ensure_ascii=False, indent=4)

    def save_experiment(self, experiment_id: ExperimentId, pdb_contents: List[str]):
        experiment_folder = self.experiment_folder(experiment_id)
        result_folder = 'results'
        for i, pdb_content in enumerate(pdb_contents):
            result_pdb_path_template = os.path.join(experiment_folder,
                                                    result_folder,
                                                    self._protein_design_file_template)
            result_pdb_path_template = result_pdb_path_template.replace('{0}', str(i))
            with open(result_pdb_path_template, 'w', encoding='utf-8') as protein_designs_f:
                protein_designs_f.write(pdb_content)

    def get_experiment_data(self, experiment_id: ExperimentId) -> List[str]:
        experiment_folder = self.experiment_folder(experiment_id)
        result_folder = 'results'
        result_pdb_path = os.path.join(experiment_folder, result_folder)

        result_pdbs = []
        for pdb_file in glob.glob(os.path.join(result_pdb_path, '*.pdb')):
            with open(pdb_file, 'r', encoding='utf-8') as protein_designs_f:
                result_pdbs.append(protein_designs_f.read())

        return result_pdbs

