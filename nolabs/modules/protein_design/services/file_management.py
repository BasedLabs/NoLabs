import glob
import json
import os.path
from io import BytesIO
from typing import List

import pathlib

from fastapi import UploadFile

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.api_models.protein_design import RunProteinDesignRequest
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.modules.file_management_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.protein_design_experiments_folder, settings.protein_design_metadata_file_name)
        self._settings = settings
        self.ensure_experiments_folder_exists()
        self._experiment_properties_filename = 'properties.json'

    def properties_exists(self, experiment_id: ExperimentId) -> bool:
        return os.path.exists(os.path.join(self.experiment_folder(experiment_id), self._experiment_properties_filename))

    async def set_properties(self, experiment_id: ExperimentId, request: RunProteinDesignRequest):
        if not request or not request.pdb_file or not request.pdb_file.filename:
            raise NoLabsException(['No PDB file provided'], error_code=ErrorCodes.protein_design_update_metadata_error)

        await request.pdb_file.seek(0)

        try:
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
        finally:
            await request.pdb_file.seek(0)

    async def set_result(self, experiment_id: ExperimentId, pdbs_content: List[str], request: RunProteinDesignRequest):
        if not request or not request.pdb_file or not request.pdb_file.filename:
            raise NoLabsException(['No PDB file provided'], error_code=ErrorCodes.protein_design_update_metadata_error)

        experiment_folder = self.experiment_folder(experiment_id)

        pdb_files = glob.glob(os.path.join(experiment_folder, '*.pdb'))
        for pdb_file in pdb_files:
            os.remove(pdb_file)

        file = pathlib.Path(request.pdb_file.filename)

        for i, pdb_content in enumerate(pdbs_content):
            protein_design_file_name = f'{file.stem}_generated_{i}{file.suffix}'
            results_pdb_path = os.path.join(experiment_folder, protein_design_file_name)
            with open(results_pdb_path, 'w', encoding='utf-8') as protein_design_f:
                protein_design_f.write(pdb_content)

        await request.pdb_file.seek(0)

    def get_properties(self, experiment_id: ExperimentId) -> RunProteinDesignRequest:
        experiment_folder = self.experiment_folder(experiment_id)

        metadata = self.get_metadata(experiment_id=experiment_id)

        with open(os.path.join(experiment_folder, self._experiment_properties_filename), 'r') as f:
            properties = json.load(f)

        return RunProteinDesignRequest(
            pdb_file=UploadFile(
                BytesIO(open(os.path.join(experiment_folder, properties['input_file_name']), 'rb').read()),
                filename=properties['input_file_name']),
            contig=properties['contig'],
            number_of_designs=properties['number_of_designs'],
            timesteps=properties['timesteps'],
            hotspots=properties['hotspots'],
            experiment_name=metadata.name.value,
            experiment_id=experiment_id.value
        )

    def get_result(self, experiment_id: ExperimentId) -> List[str]:
        experiment_folder = self.experiment_folder(experiment_id)
        result = []
        for i, pdb_file in enumerate(glob.glob(os.path.join(experiment_folder, '*_generated_*.pdb'))):
            with open(pdb_file, 'r', encoding='utf-8') as protein_design_f:
                result.append(protein_design_f.read())
        return result
