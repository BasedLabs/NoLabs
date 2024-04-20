from nolabs.api_models.biobuddy import EditMessageRequest, EditMessageResponse
from nolabs.api_models.biobuddy import Message as ApiMessage
from nolabs.api_models.biobuddy import RegularMessage as ApiRegularMessage
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.biobuddy.data_models.message import Message, RegularMessage
from nolabs.modules.biobuddy.file_management import FileManagement
from nolabs.infrastructure.settings import Settings

class EditMessageFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._file_management = file_management
        self._settings = settings

    def handle(self, request: EditMessageRequest) -> EditMessageResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        message_id = request.message_id

        self._file_management.edit_message_in_conversation(experiment_id,
                                                           message_id,
                                                            Message(id=message_id,
                                                                     role='user',
                                                                     message=RegularMessage(
                                                                         content=request.message_content
                                                                     ),
                                                                     type='text')
                                                            )

        return EditMessageResponse(ApiMessage(id=message_id,
                                              role='user',
                                              message=ApiRegularMessage(
                                                  content=request.message_content
                                              ),
                                              type='text')
                                   )
