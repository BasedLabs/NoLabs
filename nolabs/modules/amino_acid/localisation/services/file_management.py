import glob
import json
import os.path
from typing import List

from slugify import slugify

from nolabs.api_models.amino_acid.localisation import AminoAcidResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.amino_acid.file_management_base import AminoAcidFileManagementBase
from nolabs.infrastructure.settings import Settings


class FileManagement(AminoAcidFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.localisation_experiments_folder, settings.solubility_metadata_file_name)

    def set_result(self, experiment_id: ExperimentId, data: List[AminoAcidResponse]):
        experiment_folder = self.experiment_folder(experiment_id)
        for amino_acid in data:
            results_path = os.path.join(experiment_folder, f'{slugify(amino_acid.name)}_aminoacid.json')
            with open(results_path, 'w', encoding='utf-8') as localisation_f:
                localisation_f.write(json.dumps({
                    'name': amino_acid.name,
                    'sequence': amino_acid.sequence,
                    'cytosolic_proteins': amino_acid.cytosolic_proteins,
                    'mitochondrial_proteins': amino_acid.mitochondrial_proteins,
                    'nuclear_proteins': amino_acid.nuclear_proteins,
                    'other_proteins': amino_acid.other_proteins,
                    'extracellular_secreted_proteins': amino_acid.extracellular_secreted_proteins
                }))

    def get_result(self, experiment_id: ExperimentId) -> List[AminoAcidResponse]:
        experiment_folder = self.experiment_folder(experiment_id)

        results = []
        for file in glob.glob(os.path.join(experiment_folder, '*_aminoacid.json')):
            with open(file, 'r', encoding='utf-8') as amino_acid:
                j = json.loads(amino_acid.read())
                results.append(AminoAcidResponse(name=j['name'], sequence=j['sequence'],
                                                 cytosolic_proteins=j['cytosolic_proteins'],
                                                 mitochondrial_proteins=j['mitochondrial_proteins'],
                                                 nuclear_proteins=j['nuclear_proteins'],
                                                 other_proteins=j['other_proteins'],
                                                 extracellular_secreted_proteins=j['extracellular_secreted_proteins'],
                                                 )
                               )
        return results
