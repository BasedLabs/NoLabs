from typing import Dict, List

from nolabs.api_models.drug_discovery import ExperimentMetadataResponse
from nolabs.features.drug_discovery.services.file_management import FileManagement


class GetExperimentsFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self) -> List[ExperimentMetadataResponse]:
        d = self._file_management.get_all_experiments_metadata()
        result_list = []
        for experiment_metadata in self._file_management.get_all_experiments_metadata():
            result_list.append(ExperimentMetadataResponse(
                experiment_id=experiment_metadata.id.value,
                experiment_name=experiment_metadata.name.value,
                experiment_date=experiment_metadata.date
            ))

        return result_list
