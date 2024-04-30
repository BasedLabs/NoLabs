import json
import os.path
from dataclasses import asdict
from typing import Dict

from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.modules.biobuddy.data_models.message import Message, MessageType, RegularMessage, FunctionCall, \
    FunctionParam
from nolabs.modules.file_management_base import ExperimentsFileManagementBase
from nolabs.infrastructure.settings import Settings


class FileManagement(ExperimentsFileManagementBase):
    def __init__(self, settings: Settings):
        super().__init__(settings.drug_discovery_experiments_folder,
                         settings.drug_discovery_experiment_metadata_file_name)

        self._settings = settings

    def conversation_file_path(self, experiment_id: ExperimentId) -> str:
        return os.path.join(self._experiments_folder, experiment_id.value, self._settings.biobuddy_conversation_file_name)

    def ensure_conversations_file_exists(self, experiment_id: ExperimentId):
        file_path = self.conversation_file_path(experiment_id)
        if not os.path.exists(os.path.dirname(file_path)):
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
        if not os.path.exists(file_path):
            with open(file_path, 'w') as f:
                json.dump([], f)

    def append_message_to_conversation(self, experiment_id: ExperimentId, message: Message):
        self.ensure_conversations_file_exists(experiment_id)
        file_path = self.conversation_file_path(experiment_id)
        message_dict = self._message_to_dict(message=message)
        with open(file_path, 'r+') as f:
            messages = json.load(f)
            messages.append(message_dict)
            f.seek(0)
            json.dump(messages, f, indent=4)

    def edit_message_in_conversation(self, experiment_id: ExperimentId, message_id: str, new_message: Message):
        self.ensure_conversations_file_exists(experiment_id)
        file_path = self.conversation_file_path(experiment_id)
        with open(file_path, 'r+') as f:
            messages = json.load(f)
            for idx in range(len(messages)):
                message = messages[idx]
                if message['id'] == message_id:
                    message_dict = self._message_to_dict(new_message)
                    messages[idx] = message_dict
                    messages = messages[:idx+1]
                    f.seek(0)
                    json.dump(messages, f, indent=4)
                    f.truncate()
                    return

    def load_conversation(self, experiment_id: ExperimentId) -> list[Message]:
        self.ensure_conversations_file_exists(experiment_id)
        file_path = self.conversation_file_path(experiment_id)
        with open(file_path, 'r') as f:
            messages_dict = json.load(f)

        messages = []

        for message in messages_dict:
            msg = self._dict_to_message(message)
            messages.append(msg)

        return messages

    def _message_to_dict(self, message: Message) -> Dict:
        content = None
        if message.type == MessageType.FUNCTIONS:
            content = [asdict(func) for func in message.message]
        else:
            content = asdict(message.message)

        print(content)

        messsage_dict = {
            "id": message.id,
            "role": message.role,
            "message": content,
            "type": message.type.value
        }
        return messsage_dict

    def _dict_to_message(self, message: Dict) -> Message:
        if message["type"] == MessageType.TEXT.value:
            content = RegularMessage(content=message["message"]["content"])
            msg = Message(id=message["id"],
                          role=message["role"],
                          message=content,
                          type=MessageType.TEXT)
            return msg
        elif message["type"] == MessageType.FUNCTIONS.value:
            functions = []
            for func_dict in message["message"]:
                func = FunctionCall(function_name=func_dict["function_name"],
                                    parameters=[FunctionParam(name=param["name"],
                                                              value=param["value"]
                                                              ) for param in func_dict["parameters"]])
                functions.append(func)
            msg = Message(id=message["id"],
                          role=message["role"],
                          message=functions,
                          type=MessageType.FUNCTIONS)
            return msg
        else:
            raise NoLabsException(ErrorCodes.biobuddy_unexpected_message_type)