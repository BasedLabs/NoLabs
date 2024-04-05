import glob
import json
import os.path
from typing import List

from nolabs.api_models.amino_acid.solubility import AminoAcidResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.amino_acid.file_management_base import AminoAcidFileManagementBase
from nolabs.infrastructure.settings import Settings


class FileManagement(AminoAcidFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.solubility_experiments_folder, settings.solubility_metadata_file_name)

    def set_result(self, experiment_id: ExperimentId, data: List[AminoAcidResponse]):
        experiment_folder = self.experiment_folder(experiment_id)
        previous_files = glob.glob(os.path.join(experiment_folder, '*aminoacid.json'))
        for amino_acid in data:
            results_path = os.path.join(experiment_folder, f'{amino_acid.name}_aminoacid.json')
            with open(results_path, 'w', encoding='utf-8') as localisation_f:
                localisation_f.write(json.dumps({
                    'name': amino_acid.name,
                    'sequence': amino_acid.sequence,
                    'soluble_probability': amino_acid.soluble_probability
                }))
        for file in previous_files:
            os.remove(file)

    def get_result(self, experiment_id: ExperimentId) -> List[AminoAcidResponse]:
        experiment_folder = self.experiment_folder(experiment_id)

        results = []
        for file in glob.glob(os.path.join(experiment_folder, '*_aminoacid.json')):
            with open(file, 'r', encoding='utf-8') as amino_acid:
                j = json.loads(amino_acid.read())
                results.append(AminoAcidResponse(name=j['name'], sequence=j['sequence'],
                                                 soluble_probability=j['soluble_probability']))
        return results
