from nolabs.api_models.biobuddy import SaveMessageRequest, SaveMessageResponse
from nolabs.api_models.biobuddy import Message as ApiMessage
from nolabs.api_models.biobuddy import RegularMessage as ApiRegularMessage
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.biobuddy.data_models.message import Message, RegularMessage
from nolabs.modules.biobuddy.file_management import FileManagement
from nolabs.infrastructure.settings import Settings
from nolabs.utils import generate_uuid


class SendMessageFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._file_management = file_management
        self._settings = settings

    def handle(self, request: SaveMessageRequest) -> SaveMessageResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        message_id = generate_uuid()

        self._file_management.append_message_to_conversation(experiment_id,
                                                             Message(id=message_id,
                                                                     role='user',
                                                                     message=RegularMessage(
                                                                         content=request.message_content
                                                                     ),
                                                                     type='text')
                                                             )

        return SaveMessageResponse(ApiMessage(id=message_id,
                                              role='user',
                                              message=ApiRegularMessage(
                                                  content=request.message_content
                                              ),
                                              type='text')
                                   )
