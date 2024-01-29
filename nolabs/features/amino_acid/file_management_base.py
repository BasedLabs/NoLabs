import glob
import json
import os
import pathlib
from io import BytesIO

from fastapi import UploadFile

from nolabs.api_models.amino_acid.common_models import RunAminoAcidRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.features.file_management_base import ExperimentsFileManagementBase


class AminoAcidFileManagementBase(ExperimentsFileManagementBase):
    def __init__(self, experiments_folder: str, metadata_file: str):
        super().__init__(experiments_folder, metadata_file)
        self.ensure_experiments_folder_exists()
        self._experiment_properties_filename = 'properties.json'

    async def get_properties(self, experiment_id: ExperimentId) -> RunAminoAcidRequest:
        experiment_folder = self.experiment_folder(experiment_id)
        metadata = self.get_metadata(experiment_id=experiment_id)

        properties_path = os.path.join(experiment_folder, self._experiment_properties_filename)

        sequence: str | None = None

        if os.path.exists(properties_path):
            with open(properties_path, 'r') as f:
                sequence = json.load(f)['sequence']

        fastas = [
            UploadFile(
                BytesIO(open(fasta_path, 'rb').read()),
                filename=pathlib.Path(fasta_path).name
            ) for fasta_path in glob.glob(os.path.join(experiment_folder, '*.fasta'))
        ]

        return RunAminoAcidRequest(
            experiment_id=experiment_id.value,
            experiment_name=metadata.name.value,
            amino_acid_sequence=sequence,
            fastas=fastas
        )

    async def set_properties(self, experiment_id: ExperimentId, request: RunAminoAcidRequest):
        experiment_folder = self.experiment_folder(experiment_id)

        if request.amino_acid_sequence:
            properties_path = os.path.join(experiment_folder, self._experiment_properties_filename)
            with open(properties_path, 'w', encoding='utf-8') as f:
                json.dump({
                    'sequence': request.amino_acid_sequence
                }, f, ensure_ascii=False, indent=4)

        if request.fastas:
            for fasta in request.fastas:
                if not fasta.filename:
                    raise NoLabsException(['Cannot obtain name of fasta file'],
                                          ErrorCodes.amino_acid_localisation_run_error)
                fasta_content = await fasta.read()
                with open(os.path.join(experiment_folder, fasta.filename), 'wb') as f:
                    f.write(fasta_content)
                await fasta.seek(0)
