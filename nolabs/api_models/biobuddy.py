from pydantic import dataclasses as pcdataclass
from typing import List, Any, Union


@pcdataclass.dataclass
class File:
    name: str
    content: str


@pcdataclass.dataclass
class RegularMessage:
    content: str


@pcdataclass.dataclass
class FunctionParam:
    name: str
    value: Any = None


@pcdataclass.dataclass
class FunctionCall:
    function_name: str
    parameters: List[FunctionParam]


@pcdataclass.dataclass
class Message:
    role: str
    message: Union[RegularMessage,FunctionCall]
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
