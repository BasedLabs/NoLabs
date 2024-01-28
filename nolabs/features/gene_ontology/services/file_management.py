import dataclasses
import glob
import json
import os.path
from typing import List, Tuple

from nolabs.domain.amino_acid import AminoAcid
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.domain.gene_ontology import OboNode
from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings
from nolabs.utils import utcnow


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.go_experiments_folder, settings.go_metadata_file_name)
        self._settings = settings
        self._go_file_name = settings.go_file_name
        self.ensure_experiments_folder_exists()

    def set_metadata(self, experiment_id: ExperimentId, experiment_name: ExperimentName,
                     amino_acids: List[AminoAcid]):
        self.ensure_experiment_folder_exists(experiment_id)
        j = {
            'id': experiment_id.value,
            'name': experiment_name.value,
            'date': str(utcnow()),
            'properties': [dataclasses.asdict(amino_acid) for amino_acid in amino_acids]
        }

        metadata_file_path = os.path.join(self.experiment_folder(experiment_id),
                                          self._settings.solubility_metadata_file_name)
        with open(metadata_file_path, 'w', encoding='utf-8') as f:
            json.dump(j, f, ensure_ascii=False, indent=4)

    def save_experiment(self, experiment_id: ExperimentId, data: List[Tuple[AminoAcid, List[OboNode]]]):
        self.ensure_experiment_folder_exists(experiment_id)
        experiment_folder = self.experiment_folder(experiment_id)
        previous_files = glob.glob(os.path.join(experiment_folder, '*aminoacid.json'))
        for amino_acid, prob in data:
            results_path = os.path.join(experiment_folder, f'{amino_acid.name}_aminoacid.json')
            with open(results_path, 'w', encoding='utf-8') as go_f:
                go_f.write(json.dumps({
                    'name': amino_acid.name,
                    'sequence': amino_acid.sequence,
                    'go': [
                        dataclasses.asdict(g) for g in prob
                    ]
                }))
        for file in previous_files:
            os.remove(file)

    def get_experiment_data(self, experiment_id: ExperimentId) -> List[Tuple[AminoAcid, List[OboNode]]]:
        experiment_folder = self.experiment_folder(experiment_id)

        results = []
        for file in glob.glob(os.path.join(experiment_folder, '*_aminoacid.json')):
            with open(file, 'r', encoding='utf-8') as amino_acid:
                j = json.loads(amino_acid.read())
                results.append(
                    (AminoAcid(name=j['name'],
                               sequence=j['sequence']),
                     [OboNode(
                         name=go['name'],
                         namespace=go['namespace'],
                         edges=go['edges']
                     ) for go in j['go']]
                     ))
        return results
