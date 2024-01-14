import dataclasses
import glob
import json
import os.path
from typing import List, Tuple

from nolabs.domain.amino_acid import AminoAcid
from nolabs.features.file_magement_base import ExperimentsFileManagementBase
from nolabs.domain.localisation import LocalisationProbability
from nolabs.domain.experiment import ExperimentId, ExperimentName, ExperimentMetadata
from nolabs.infrastructure.settings import Settings
from nolabs.utils.datetime_utils import utcnow


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.solubility_experiments_folder, settings.solubility_metadata_file_name)
        self._settings = settings
        self.ensure_experiments_folder_exists()

    def update_metadata(self, experiment_id: ExperimentId, experiment_name: ExperimentName,
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

    def save_experiment(self, experiment_id: ExperimentId, data: List[Tuple[AminoAcid, LocalisationProbability]]):
        self.ensure_experiment_folder_exists(experiment_id)
        experiment_folder = self.experiment_folder(experiment_id)
        previous_files = glob.glob(os.path.join(experiment_folder, '*aminoacid.json'))
        for amino_acid, prob in data:
            results_path = os.path.join(experiment_folder, f'{amino_acid.name}_aminoacid.json')
            with open(results_path, 'w', encoding='utf-8') as localisation_f:
                localisation_f.write(json.dumps({
                    'name': amino_acid.name,
                    'sequence': amino_acid.sequence,
                    **dataclasses.asdict(prob)
                }))
        for file in previous_files:
            os.remove(file)

    def get_experiment_data(self, experiment_id: ExperimentId) -> List[Tuple[AminoAcid, LocalisationProbability]]:
        experiment_folder = self.experiment_folder(experiment_id)

        results = []
        for file in glob.glob(os.path.join(experiment_folder, '*_aminoacid.json')):
            with open(file, 'r', encoding='utf-8') as amino_acid:
                j = json.loads(amino_acid.read())
                results.append((AminoAcid(name=j['name'], sequence=j['sequence']),
                                LocalisationProbability(
                                    nuclear_proteins=j['nuclear_proteins'],
                                    other_proteins=j['other_proteins'],
                                    mitochondrial_proteins=j['mitochondrial_proteins'],
                                    cytosolic_proteins=j['cytosolic_proteins'],
                                    extracellular_secreted_proteins=j['extracellular_secreted_proteins']
                                )
                                ))
        return results
