import glob
import json
import os.path
import pathlib
from io import BytesIO
from typing import List, Tuple

from dateutil import parser
from fastapi import UploadFile

from nolabs.api_models.conformations import RunSimulationsRequest, TimelineResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.modules.file_management_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.conformations_experiments_folder, settings.conformations_metadata_file_name)
        self._experiment_properties_filename = 'properties.json'
        self._timeline_path = 'timeline.json'

    def properties_exists(self, experiment_id: ExperimentId) -> bool:
        return os.path.exists(os.path.join(self.experiment_folder(experiment_id), self._experiment_properties_filename))

    async def set_properties(self, experiment_id: ExperimentId, request: RunSimulationsRequest):
        if not request or not request.pdb_file or not request.pdb_file.filename:
            raise NoLabsException(['No PDB file provided'], error_code=ErrorCodes.conformations_update_metadata_error)

        self.ensure_experiment_folder_exists(experiment_id=experiment_id)

        await request.pdb_file.seek(0)

        try:
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
        finally:
            await request.pdb_file.seek(0)

    async def set_result(self,
                         experiment_id: ExperimentId,
                         pdb_content: str | None,
                         timeline: List[TimelineResponse],
                         request: RunSimulationsRequest):
        if not request or not request.pdb_file or not request.pdb_file.filename:
            raise NoLabsException(['No PDB file provided'], error_code=ErrorCodes.conformations_update_metadata_error)

        self.ensure_experiment_folder_exists(experiment_id=experiment_id)

        await request.pdb_file.seek(0)

        try:
            experiment_folder = self.experiment_folder(experiment_id)

            if pdb_content:
                file = pathlib.Path(request.pdb_file.filename)

                protein_design_file_name = f'{file.stem}_generated{file.suffix}'
                results_pdb_path = os.path.join(experiment_folder, protein_design_file_name)
                with open(results_pdb_path, 'w', encoding='utf-8') as protein_design_f:
                    protein_design_f.write(pdb_content)

            timeline_j = []
            for tm in timeline:
                timeline_j.append({
                    'message': tm.message,
                    'error': tm.error,
                    'created_at': str(tm.created_at)
                })

            timeline_path = os.path.join(experiment_folder, self._timeline_path)
            with open(timeline_path, 'w', encoding='utf-8') as timeline_f:
                json.dump(timeline_j, timeline_f)
        finally:
            await request.pdb_file.seek(0)

    def get_properties(self, experiment_id: ExperimentId) -> RunSimulationsRequest:
        experiment_folder = self.experiment_folder(experiment_id)

        with open(os.path.join(experiment_folder, self._experiment_properties_filename), 'r') as f:
            properties = json.load(f)

        metadata = self.get_metadata(experiment_id=experiment_id)

        return RunSimulationsRequest(
            pdb_file=UploadFile(
                BytesIO(open(os.path.join(experiment_folder, properties['input_file_name']), 'rb').read()),
                filename=properties['input_file_name']),
            total_frames=properties['total_frames'],
            temperature_k=properties['temperature_k'],
            take_frame_every=properties['take_frame_every'],
            step_size=properties['step_size'],
            replace_non_standard_residues=properties['replace_non_standard_residues'],
            add_missing_atoms=properties['add_missing_atoms'],
            add_missing_hydrogens=properties['add_missing_hydrogens'],
            friction_coeff=properties['friction_coeff'],
            ignore_missing_atoms=properties['ignore_missing_atoms'],
            integrator=properties['integrator'],
            experiment_id=experiment_id.value,
            experiment_name=metadata.name.value
        )

    def get_result(self, experiment_id: ExperimentId) -> Tuple[str | None, List[TimelineResponse]]:
        experiment_folder = self.experiment_folder(experiment_id)
        protein_content = None
        for i, pdb_file in enumerate(glob.glob(os.path.join(experiment_folder, '*_generated*.pdb'))):
            with open(pdb_file, 'r', encoding='utf-8') as conformations_f:
                protein_content = conformations_f.read()

        timeline_path = os.path.join(experiment_folder, self._timeline_path)
        if os.path.exists(timeline_path):
            timeline = json.load(open(timeline_path, 'r'))
        else:
            timeline = []

        return (protein_content, [TimelineResponse(
            message=tl['message'],
            error=tl['error'],
            created_at=parser.parse(tl['created_at'])
        )
            for tl in timeline])
