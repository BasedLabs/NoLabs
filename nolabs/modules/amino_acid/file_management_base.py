import glob
import json
import os
import pathlib
from io import BytesIO

import slugify
from fastapi import UploadFile

from nolabs.api_models.amino_acid.common_models import RunAminoAcidRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.modules.file_management_base import ExperimentsFileManagementBase


class AminoAcidFileManagementBase(ExperimentsFileManagementBase):
    def __init__(self, experiments_folder: str, metadata_file: str):
        super().__init__(experiments_folder, metadata_file)
        self.ensure_experiments_folder_exists()
        self._experiment_properties_filename = 'properties.json'

    def properties_exists(self, experiment_id: ExperimentId) -> bool:
        return os.path.exists(os.path.join(self.experiment_folder(experiment_id), self._experiment_properties_filename))

    async def get_properties(self, experiment_id: ExperimentId) -> RunAminoAcidRequest:
        experiment_folder = self.experiment_folder(experiment_id)
        metadata = self.get_metadata(experiment_id=experiment_id)

        properties_path = os.path.join(experiment_folder, self._experiment_properties_filename)

        with open(properties_path, 'r') as f:
            sequences = json.load(f)

        amino_acid = [s for s in sequences if s['filename'] == 'unknown.fasta']
        fastas = [s for s in sequences if s['filename'] != 'unknown.fasta']

        return RunAminoAcidRequest(
            experiment_id=experiment_id.value,
            experiment_name=metadata.name.value,
            amino_acid_sequence=amino_acid[0] if amino_acid else None,
            fastas=[
                UploadFile(
                    BytesIO(bytes(f['sequence'], 'utf-8')),
                    filename=f['filename']
                )
                for f in fastas]
        )

    async def set_properties(self, experiment_id: ExperimentId, request: RunAminoAcidRequest):
        sequences = []

        if request.amino_acid_sequence:
            sequences.append({'sequence': request.amino_acid_sequence, 'filename': 'unknown.fasta'})

        if request.fastas:
            for fasta in request.fastas:
                if not fasta.filename:
                    raise NoLabsException(['Cannot obtain name of fasta file'],
                                          ErrorCodes.amino_acid_localisation_run_error)
                sequences.append(
                    {'sequence': (await fasta.read()).decode('utf-8'), 'filename': slugify.slugify(fasta.filename)})
                await fasta.seek(0)

        experiment_folder = self.experiment_folder(experiment_id)
        properties_path = os.path.join(experiment_folder, self._experiment_properties_filename)
        with open(properties_path, 'w', encoding='utf-8') as f:
            json.dump(sequences, f, ensure_ascii=False, indent=4)
