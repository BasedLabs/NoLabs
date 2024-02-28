from nolabs.api_models.biobuddy import SendMessageRequest, SendMessageResponse
from nolabs.api_models.biobuddy import Message as ApiMessage
from nolabs.domain.experiment import ExperimentId
from nolabs.features.biobuddy.data_models.message import Message
from nolabs.features.biobuddy.file_management import FileManagement

import biobuddy_microservice
from biobuddy_microservice import DefaultApi
from biobuddy_microservice import ApiClient
from biobuddy_microservice import Configuration

from nolabs.infrastructure.settings import Settings


class SendMessageFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._file_management = file_management
        self._settings = settings

    def handle(self, request: SendMessageRequest) -> SendMessageResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        #messages = self._file_management.load_conversation(experiment_id)

        configuration = Configuration(
            host=self._settings.biobuddy_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            request = biobuddy_microservice.SendMessageToBioBuddyRequest(message_content=request.message_content, previous_messages=[])
            reply = api_instance.predict_send_message_post(send_message_to_bio_buddy_request=request).chatgpt_reply
            self._file_management.update_conversation(experiment_id,
                                                      Message(role='user', content=request.message_content, type='text'))
            self._file_management.update_conversation(experiment_id, Message(role='assistant', content=reply, type='text'))

        return SendMessageResponse(biobuddy_response=ApiMessage(role='user', content=reply, type='text'))
