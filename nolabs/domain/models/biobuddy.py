__all__ = ["Chat", "Message", "TextMessage", "FunctionCall", "UserRoleEnum"]

from enum import Enum
from typing import Any, List, Union
from uuid import UUID

from mongoengine import (Document, EmbeddedDocument, EmbeddedDocumentListField,
                         EnumField, StringField, UUIDField)


class UserRoleEnum(str, Enum):
    user = "user"
    biobuddy = "biobuddy"


class Message(EmbeddedDocument):
    """
    Base class for all messages
    """

    message_id: UUID = UUIDField(required=True)
    sender: UserRoleEnum = EnumField(UserRoleEnum, required=True)

    meta = {"allow_inheritance": True}


class TextMessage(Message):
    """
    A regular text message
    """

    content: str = StringField(required=True)


class FunctionParam(EmbeddedDocument):
    name: str = StringField(required=True)
    value: Any = StringField(required=False)


class FunctionCall(EmbeddedDocument):
    function_name: str = StringField(required=True)
    arguments = EmbeddedDocumentListField(FunctionParam)


class FunctionCallMessage(Message):
    """
    A message that contains function calls
    """

    function_calls: List[FunctionCall] = EmbeddedDocumentListField(FunctionCall)


class Chat(Document):
    """
    Chat document containing messages between user and biobuddy
    """

    experiment_id: UUID = UUIDField(required=True)
    messages: List[Union[TextMessage, FunctionCallMessage]] = EmbeddedDocumentListField(
        Message
    )

    def add_text_message(self, message_id: UUID, sender: UserRoleEnum, content: str):
        if sender not in [UserRoleEnum.user, UserRoleEnum.biobuddy]:
            raise ValueError("Invalid sender")

        if isinstance(content, str):
            message = TextMessage(message_id=message_id, sender=sender, content=content)
            self.messages.append(message)
        else:
            raise ValueError("Invalid message content for text message")

    def load_conversation(self, experiment_id: UUID) -> "Chat":
        """
        Load the conversation for the given experiment_id
        """
        return Chat.objects(experiment_id=experiment_id).first()

    def edit_message(self, message_id: UUID, new_content: str):
        """
        Edit a user message. All messages after the edited one will be removed.
        """
        for i, message in enumerate(self.messages):
            if message.message_id == message_id:
                if message.sender != UserRoleEnum.user:
                    raise ValueError("Only user messages can be edited")
                if isinstance(message, TextMessage):
                    message.content = new_content
                    self.messages = self.messages[: i + 1]
                    return
                else:
                    raise ValueError("Message content type mismatch")
        raise ValueError("Message ID not found")

    def add_function_call_message(
        self, message_id: UUID, sender: UserRoleEnum, function_calls: List[FunctionCall]
    ):
        if sender != UserRoleEnum.biobuddy:
            raise ValueError("Only biobuddy can send function call messages")

        if isinstance(function_calls, list):
            message = FunctionCallMessage(
                message_id=message_id, sender=sender, function_calls=function_calls
            )
            self.messages.append(message)
        else:
            raise ValueError("Invalid message content for function call message")
