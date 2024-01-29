import dataclasses
import glob
import json
import os.path
import pathlib
from io import BytesIO
from typing import List, Tuple

from fastapi import UploadFile
from slugify import slugify

from nolabs.api_models.gene_ontology import AminoAcidResponse, RunGeneOntologyResponseDataNode, RunGeneOntologyRequest
from nolabs.domain.amino_acid import AminoAcid
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.domain.gene_ontology import OboNode
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings
from nolabs.utils import utcnow


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.go_experiments_folder, settings.go_metadata_file_name)
        self._settings = settings
        self._go_file_name = settings.go_file_name
        self.ensure_experiments_folder_exists()
        self._experiment_properties_filename = 'properties.json'

    def set_result(self, experiment_id: ExperimentId, data: List[AminoAcidResponse]):
        self.ensure_experiment_folder_exists(experiment_id)
        experiment_folder = self.experiment_folder(experiment_id)
        for amino_acid in data:
            results_path = os.path.join(experiment_folder, f'{slugify(amino_acid.name)}_aminoacid.json')
            with open(results_path, 'w', encoding='utf-8') as go_f:
                go_f.write(json.dumps({
                    'name': amino_acid.name,
                    'sequence': amino_acid.sequence,
                    'go': {
                        go_name: {
                            'name': node.name,
                            'namespace': node.namespace,
                            'edges': node.edges
                        }
                        for go_name, node in amino_acid.go.items()
                    }
                }))

    async def set_properties(self, experiment_id: ExperimentId, request: RunGeneOntologyRequest):
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

    async def get_properties(self, experiment_id: ExperimentId) -> RunGeneOntologyRequest:
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

        return RunGeneOntologyRequest(
            experiment_id=experiment_id.value,
            experiment_name=metadata.name.value,
            amino_acid_sequence=sequence,
            fastas=fastas
        )

    def get_result(self, experiment_id: ExperimentId) -> List[AminoAcidResponse]:
        experiment_folder = self.experiment_folder(experiment_id)

        results = []
        for file in glob.glob(os.path.join(experiment_folder, '*_aminoacid.json')):
            with open(file, 'r', encoding='utf-8') as amino_acid:
                j = json.loads(amino_acid.read())
                results.append(
                    AminoAcidResponse(name=j['name'],
                                      sequence=j['sequence'],
                                      go={go_name: RunGeneOntologyResponseDataNode(
                                          name=node['name'],
                                          namespace=node['namespace'],
                                          edges=node['edges']
                                      ) for go_name, node in j['go'].items()}
                                      )
                )
        return results
