from nolabs.api_models.experiment import ExperimentMetadataResponse
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.features.file_management_base import ExperimentsFileManagementBase
from nolabs.utils import uuid_utils, datetime_utils


class CreateExperimentFeature:
    def __init__(self, file_management: ExperimentsFileManagementBase):
        self._file_management = file_management

    async def handle(self) -> ExperimentMetadataResponse:
        await self._file_management.set_metadata(experiment_id=ExperimentId(uuid_utils.generate_uuid()),
                                           experiment_name=ExperimentName('New experiment'))

        return ExperimentMetadataResponse(
            id=uuid_utils.generate_uuid(),
            name='New experiment',
            date=datetime_utils.utcnow()
        )
