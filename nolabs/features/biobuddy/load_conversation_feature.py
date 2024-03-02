from nolabs.api_models.biobuddy import LoadConversationRequest, LoadConversationResponse
from nolabs.api_models.biobuddy import Message as ApiMessage
from nolabs.domain.experiment import ExperimentId
from nolabs.features.biobuddy.file_management import FileManagement


class LoadConversationFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, request: LoadConversationRequest) -> LoadConversationResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        messages = self._file_management.load_conversation(experiment_id)

        return LoadConversationResponse(messages=[ApiMessage(role=message.role,
                                                             content=message.content,
                                                             type=message.type) for message in messages])
