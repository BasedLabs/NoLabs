import glob
import json
import os.path
from typing import List

import pathlib

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.api_models.protein_design import RunProteinDesignRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.domain.protein_design import ExperimentProperties
from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings
from nolabs.utils import utcnow


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__('protein_design', 'metadata.json')
        self._settings = settings
        self.ensure_experiments_folder_exists()
        self._experiment_properties_filename = 'properties.json'

    async def save_experiment(self, experiment_id: ExperimentId, pdbs_content: List[str], request: RunProteinDesignRequest):
        if not request or not request.pdb_file or not request.pdb_file.filename:
            raise NoLabsException(['No PDB file provided'], error_code=ErrorCodes.protein_design_update_metadata_error)

        properties = {
            'contig': request.contig,
            'number_of_designs': request.number_of_designs,
            'input_file_name': request.pdb_file.filename,
            'timesteps': request.timesteps,
            'hotspots': request.hotspots
        }

        experiment_folder = self.experiment_folder(experiment_id)

        pdb_files = glob.glob(os.path.join(experiment_folder, '*.pdb'))
        for pdb_file in pdb_files:
            os.remove(pdb_file)

        with open(os.path.join(experiment_folder, self._experiment_properties_filename), 'w') as f:
            json.dump(properties, f, ensure_ascii=False, indent=4)

        with open(os.path.join(experiment_folder, request.pdb_file.filename), 'wb') as f:
            f.write(await request.pdb_file.read())

        file = pathlib.Path(request.pdb_file.filename)

        for i, pdb_content in enumerate(pdbs_content):
            protein_design_file_name = f'{file.stem}_generated_{i}{file.suffix}'
            results_pdb_path = os.path.join(experiment_folder, protein_design_file_name)
            with open(results_pdb_path, 'w', encoding='utf-8') as protein_design_f:
                protein_design_f.write(pdb_content)

    def experiment_properties(self, experiment_id) -> ExperimentProperties:
        experiment_folder = self.experiment_folder(experiment_id)
        with open(os.path.join(experiment_folder, self._experiment_properties_filename), 'r') as f:
            properties = json.load(f)

        with open(os.path.join(experiment_folder, properties['input_file_name']), 'r') as f:
            input_pdb_file = f.read()

        return ExperimentProperties(
            contig=properties['contig'],
            number_of_designs=properties['number_of_designs'],
            input_file_name=properties['input_file_name'],
            timesteps=properties['timesteps'],
            hotspots=properties['hotspots'],
            input_file_content=input_pdb_file
        )

    def get_experiment_data(self, experiment_id: ExperimentId) -> List[str]:
        experiment_folder = self.experiment_folder(experiment_id)
        result = []
        for i, pdb_file in enumerate(glob.glob(os.path.join(experiment_folder, '*_generated_*.pdb'))):
            with open(pdb_file, 'r', encoding='utf-8') as protein_design_f:
                result.append(protein_design_f.read())
        return result
