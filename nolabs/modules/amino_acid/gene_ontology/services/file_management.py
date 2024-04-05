import glob
import json
import os.path
from typing import List

from slugify import slugify

from nolabs.api_models.amino_acid.gene_ontology import AminoAcidResponse, RunGeneOntologyResponseDataNode
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.amino_acid.file_management_base import AminoAcidFileManagementBase
from nolabs.infrastructure.settings import Settings


class FileManagement(AminoAcidFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.go_experiments_folder, settings.go_metadata_file_name)

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
