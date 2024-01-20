import glob
import json
import os.path
from typing import List

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.api_models.protein_design import RunProteinDesignRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings
from nolabs.utils import utcnow


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.protein_design_experiments_folder, settings.protein_design_metadata_file_name)
        self._settings = settings
        self.ensure_experiments_folder_exists()
        self._experiment_properties_filename = 'properties.json'

    async def update_metadata(self, experiment_id: ExperimentId, experiment_name: ExperimentName,
                              run_protein_design_request: RunProteinDesignRequest):
        if not run_protein_design_request.pdb_file or not run_protein_design_request.pdb_file.filename:
            raise NoLabsException('No PDB file', ErrorCodes.protein_design_update_metadata_error)

        self.ensure_experiment_folder_exists(experiment_id)
        j = {
            'id': experiment_id.value,
            'name': experiment_name.value,
            'date': str(utcnow())
        }

        with open(os.path.join(self.experiment_folder(experiment_id), run_protein_design_request.pdb_file.filename),
                  'wb') as f:
            f.write(await run_protein_design_request.pdb_file.read())

        metadata_file_path = os.path.join(self.experiment_folder(experiment_id),
                                          self._settings.conformations_metadata_file_name)
        with open(metadata_file_path, 'w', encoding='utf-8') as f:
            json.dump(j, f, ensure_ascii=False, indent=4)

        properties = {
            'filename': run_protein_design_request.pdb_file.filename,
            'contig': run_protein_design_request.contig,
            'hotspots': run_protein_design_request.hotspots,
            ''
        }

    def save_experiment(self, experiment_id: ExperimentId, pdbs_content: List[str]):
        metadata = self.get_experiment_metadata(experiment_id)

        experiment_folder = self.experiment_folder(experiment_id)
        for i, pdb_content in enumerate(pdbs_content):
            protein_design_file_name = self._protein_design_file_name.replace('{0}', str(i))
            results_pdb_path = os.path.join(experiment_folder, protein_design_file_name)
            with open(results_pdb_path, 'w', encoding='utf-8') as protein_design_f:
                protein_design_f.write(pdb_content)

    def get_experiment_data(self, experiment_id: ExperimentId) -> List[str]:
        experiment_folder = self.experiment_folder(experiment_id)
        result = []
        for i, pdb_file in enumerate(glob.glob(os.path.join(experiment_folder, '*.pdb'))):
            with open(pdb_file, 'r', encoding='utf-8') as protein_design_f:
                result.append(protein_design_f.read())
        return result
