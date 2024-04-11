import json
import os.path
from dataclasses import asdict

from nolabs.domain.experiment import ExperimentId
from nolabs.modules.biobuddy.data_models.message import Message
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

    def update_conversation(self, experiment_id: ExperimentId, message: Message):
        self.ensure_conversations_file_exists(experiment_id)
        file_path = self.conversation_file_path(experiment_id)
        with open(file_path, 'r+') as f:
            messages = json.load(f)
            messages.append(asdict(message))
            f.seek(0)
            json.dump(messages, f, indent=4)

    def load_conversation(self, experiment_id: ExperimentId) -> list[Message]:
        self.ensure_conversations_file_exists(experiment_id)
        file_path = self.conversation_file_path(experiment_id)
        with open(file_path, 'r') as f:
            messages_dict = json.load(f)
        return [Message(**msg) for msg in messages_dict]
