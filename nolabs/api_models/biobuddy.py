from pydantic import dataclasses as pcdataclass
from fastapi import UploadFile
from typing import List, Optional

@pcdataclass.dataclass
class File:
    name: str
    content: str

@pcdataclass.dataclass
class Message:
    role: str
    content: str | File
    type: str

@pcdataclass.dataclass
class FunctionCallResponse:
    function_name: str
    parameters: str

@pcdataclass.dataclass
class GetExperimentRequest:
    experiment_id: str

@pcdataclass.dataclass
class LoadConversationRequest:
    experiment_id: str

@pcdataclass.dataclass
class LoadConversationResponse:
    messages: List[Message]

@pcdataclass.dataclass
class SendMessageRequest:
    experiment_id: str
    message_content: str

@pcdataclass.dataclass
class SendMessageResponse:
    biobuddy_response: Message | FunctionCallResponse

