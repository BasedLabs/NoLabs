import glob
import json
import os.path
import pathlib

from nolabs.api_models.conformations import RunSimulationsRequest
from nolabs.domain.conformations import ExperimentProperties
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__('protein_design', 'metadata.json')
        self._settings = settings
        self.ensure_experiments_folder_exists()
        self._experiment_properties_filename = 'properties.json'

    async def save_experiment(self, experiment_id: ExperimentId, pdb_content: str,
                              request: RunSimulationsRequest):
        if not request or not request.pdb_file or not request.pdb_file.filename:
            raise NoLabsException(['No PDB file provided'], error_code=ErrorCodes.conformations_update_metadata_error)

        properties = {
            'input_file_name': request.pdb_file.filename,
            'total_frames': request.total_frames,
            'temperature_k': request.temperature_k,
            'take_frame_every': request.take_frame_every,
            'step_size': request.step_size,
            'replace_non_standard_residues': request.replace_non_standard_residues,
            'add_missing_atoms': request.add_missing_atoms,
            'add_missing_hydrogens': request.add_missing_hydrogens,
            'friction_coeff': request.friction_coeff,
            'ignore_missing_atoms': request.ignore_missing_atoms,
            'integrator': request.integrator.value
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

        protein_design_file_name = f'{file.stem}_generated{file.suffix}'
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
            input_file_name=properties['input_file_name'],
            input_file=input_pdb_file,
            total_frames=properties['total_frames'],
            temperature_k=properties['temperature_k'],
            take_frame_every=properties['take_frame_every'],
            step_size=properties['step_size'],
            replace_non_standard_residues=properties['replace_non_standard_residues'],
            add_missing_atoms=properties['add_missing_atoms'],
            add_missing_hydrogens=properties['add_missing_hydrogens'],
            friction_coeff=properties['friction_coeff'],
            ignore_missing_atoms=properties['ignore_missing_atoms'],
            integrator=properties['integrator']
        )

    def get_experiment_data(self, experiment_id: ExperimentId) -> str:
        experiment_folder = self.experiment_folder(experiment_id)
        result = []
        for i, pdb_file in enumerate(glob.glob(os.path.join(experiment_folder, '*_generated*.pdb'))):
            with open(pdb_file, 'r', encoding='utf-8') as protein_design_f:
                result.append(protein_design_f.read())
        return result[0]
