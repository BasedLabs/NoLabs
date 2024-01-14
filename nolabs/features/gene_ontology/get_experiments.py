from typing import Dict, List

from nolabs.api_models.gene_ontology import ExperimentMetadataResponse
from nolabs.features.gene_ontology.services.file_management import FileManagement


class GetExperimentsFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self) -> List[ExperimentMetadataResponse]:
        d = self._file_management.get_all_experiments_metadata()
        result_d = []
        for metadata in d:
            result_d.append(ExperimentMetadataResponse(
                id=metadata.id.value,
                name=metadata.name.value,
                date=metadata.date
            ))

        return result_d
