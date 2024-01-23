from nolabs.api_models.experiment import ExperimentMetadataResponse
from nolabs.utils import uuid_utils, datetime_utils


class CreateExperimentFeature:
    def handle(self) -> ExperimentMetadataResponse:
        return ExperimentMetadataResponse(
            id=uuid_utils.generate_uuid(),
            name='New experiment',
            date=datetime_utils.utcnow()
        )
