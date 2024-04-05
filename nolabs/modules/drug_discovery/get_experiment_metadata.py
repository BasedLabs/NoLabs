from nolabs.api_models.experiment import ExperimentMetadataRequest, ExperimentMetadataResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.file_management_base import ExperimentsFileManagementBase


class GetExperimentMetaDataFeature:
    def __init__(self, file_management: ExperimentsFileManagementBase):
        self._file_management = file_management

    def handle(self, request: ExperimentMetadataRequest) -> ExperimentMetadataResponse:
        assert request

        experiment_id = ExperimentId(request.id)

        metadata = self._file_management.get_metadata(experiment_id)

        return ExperimentMetadataResponse(id=metadata.id.value, name=metadata.name.value, date=metadata.date)

