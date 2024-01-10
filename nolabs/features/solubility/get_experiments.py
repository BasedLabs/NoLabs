from typing import Dict

from nolabs.api_models.solubility import ExperimentMetadataResponse
from nolabs.features.solubility.services.file_management import FileManagement


class GetExperimentsFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self) -> Dict[str, ExperimentMetadataResponse]:
        d = self._file_management.get_all_experiments_metadata()
        result_d = {}
        for experiment_id in d.keys():
            metadata = d[experiment_id]
            result_d[experiment_id.value] = ExperimentMetadataResponse(
                id=metadata.id.value,
                name=metadata.name.value,
                date=metadata.date
            )

        return result_d
