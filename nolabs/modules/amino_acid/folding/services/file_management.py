import glob
import json
import os.path
from typing import List

from slugify import slugify

from nolabs.api_models.amino_acid.folding import AminoAcidResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.amino_acid.file_management_base import AminoAcidFileManagementBase
from nolabs.infrastructure.settings import Settings


class FileManagement(AminoAcidFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.folding_experiments_folder, settings.folding_metadata_file_name)

    def set_result(self, experiment_id: ExperimentId, data: List[AminoAcidResponse]):
        experiment_folder = self.experiment_folder(experiment_id)
        previous_files = glob.glob(os.path.join(experiment_folder, '*aminoacid.json'))
        for amino_acid in data:
            results_path = os.path.join(experiment_folder, f'{slugify(amino_acid.name)}_aminoacid.json')
            pdb_file = os.path.join(experiment_folder, amino_acid.pdb_file_name)
            with open(pdb_file, 'w') as f:
                f.write(amino_acid.pdb_file)
            with open(results_path, 'w', encoding='utf-8') as localisation_f:
                localisation_f.write(json.dumps({
                    'name': amino_acid.name,
                    'sequence': amino_acid.sequence,
                    'file_name': amino_acid.pdb_file_name
                }))
        for file in previous_files:
            os.remove(file)

    def get_result(self, experiment_id: ExperimentId) -> List[AminoAcidResponse]:
        experiment_folder = self.experiment_folder(experiment_id)

        results = []
        for file in glob.glob(os.path.join(experiment_folder, '*_aminoacid.json')):
            with open(file, 'r', encoding='utf-8') as amino_acid:
                j = json.loads(amino_acid.read())
                with open(os.path.join(experiment_folder, j['file_name']), 'r', encoding='utf-8') as f:
                    pdb_content = f.read()
                results.append(AminoAcidResponse(name=j['name'],
                                                 sequence=j['sequence'],
                                                 pdb_file=pdb_content,
                                                 pdb_file_name=j['file_name']
                                                 )
                               )
        return results
