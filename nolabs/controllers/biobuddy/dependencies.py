from typing import Annotated

from fastapi import Depends

from nolabs.controllers.common_dependencies import settings_dependency
from nolabs.features.biobuddy.file_management import FileManagement
from nolabs.features.biobuddy.load_conversation_feature import LoadConversationFeature
from nolabs.features.biobuddy.send_message_feature import SendMessageFeature
from nolabs.infrastructure.settings import Settings


def file_management_dependency(settings: Annotated[Settings, Depends(settings_dependency)]) -> FileManagement:
    return FileManagement(settings=settings)


def load_conversation_dependency(file_management: Annotated[FileManagement,
Depends(file_management_dependency)]) -> LoadConversationFeature:
    return LoadConversationFeature(file_management=file_management)


def send_message_dependency(settings: Annotated[Settings, Depends(settings_dependency)],
                            file_management: Annotated[FileManagement,
                            Depends(file_management_dependency)]) -> SendMessageFeature:
    return SendMessageFeature(settings=settings, file_management=file_management)
