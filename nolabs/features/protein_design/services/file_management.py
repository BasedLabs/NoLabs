import dataclasses
import glob
import json
import os.path
import shutil
from typing import Dict, List

from nolabs.api_models.protein_design import RunProteinDesignRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName, ExperimentMetadata
from nolabs.infrastructure.settings import Settings
from nolabs.utils.datetime_utils import DateTimeUtils


class FileManagement:
    def __init__(self, settings: Settings, dt_utils: DateTimeUtils):
        self._settings = settings
        self._dt_utils = dt_utils
        self._protein_design_file_template = settings.protein_design_file_name
        self._results_folder = 'results'

    def ensure_experiments_folder_exists(self):
        if not os.path.isdir(self._settings.protein_design_experiments_folder):
            os.mkdir(self._settings.protein_design_experiments_folder)

    def experiment_folder(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._settings.protein_design_experiments_folder, experiment_id.value)

    def create_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        if not os.path.isdir(experiment_folder):
            os.mkdir(experiment_folder)

    def delete_experiment_folder(self, experiment_id: ExperimentId):
        experiment_folder = self.experiment_folder(experiment_id)
        shutil.rmtree(experiment_folder, ignore_errors=True)

    def get_experiment_metadata(self) -> ExperimentMetadata:
        metadata_file = os.path.join(self._settings.protein_design_experiments_folder,
                                     self._settings.protein_design_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        id = ExperimentId(metadata['id'])
        return ExperimentMetadata(
            id=id,
            name=ExperimentName(metadata['name']),
            date=metadata['date'],
            properties=metadata['properties']
        )

    def get_all_experiments_metadata(self) -> Dict[ExperimentId, ExperimentMetadata]:
        metadata_files = os.path.join(self._settings.protein_design_experiments_folder,
                                      f"*{self._settings.protein_design_metadata_file_name}")
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
                        run_protein_designs_request: RunProteinDesignRequest):
        j = {
            'id': experiment_id,
            'name': experiment_name,
            'date': self._dt_utils.utcnow(),
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

    def change_experiment_name(self, experiment_id: ExperimentId, experiment_name: ExperimentName):
        metadata_file = os.path.join(self._settings.protein_design_experiments_folder,
                                     experiment_id.value,
                                     self._settings.protein_design_metadata_file_name)
        metadata = json.load(open(metadata_file, 'r', encoding='utf-8'))
        metadata['name'] = experiment_name.value
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, ensure_ascii=False, indent=4)

