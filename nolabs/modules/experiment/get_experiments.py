from typing import List

from nolabs.modules.file_management_base import ExperimentsFileManagementBase
from nolabs.api_models.experiment import ExperimentMetadataResponse


class GetExperimentsFeature:
    def __init__(self, file_management: ExperimentsFileManagementBase):
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
