from typing import Any, Dict, List, Optional, Union
from uuid import UUID

from pydantic import BaseModel
from pydantic.dataclasses import dataclass


@dataclass
class CheckBioBuddyEnabledResponse:
    enabled: bool


@dataclass
class File:
    name: str
    content: str


@dataclass
class RegularMessage:
    content: str


@dataclass
class FunctionParam:
    name: str
    value: Any = None


class FileData(BaseModel):
    content: Optional[str] = None  # Assuming file content is a string; adjust as needed
    metadata: Dict[str, Any]  # Metadata related to the file


# Subclass for RcsbPdbData
@dataclass
class RcsbPdbMetaData:
    link: str


class RcsbPdbData(FileData):
    metadata: RcsbPdbMetaData  # Example of a RcsbPdbMetaData


# Subclass for ChemBLData
@dataclass
class ChemBLMetaData:
    chembl_id: str
    link: str
    pref_name: str


class ChemBLData(FileData):
    metadata: ChemBLMetaData  # Example of a specific attribute for ChemBLData


# Base class for return data, which now contains a list of FileData objects
class FunctionCallReturnData(BaseModel):
    files: List[
        Union[RcsbPdbData, ChemBLData]
    ]  # List of files with their content and metadata


@dataclass
class FunctionCall:
    function_name: str
    arguments: List[FunctionParam]
    data: FunctionCallReturnData = None


@dataclass
class Message:
    id: UUID
    role: str
    message: Union[RegularMessage, List[FunctionCall]]
    type: str


@dataclass
class GetExperimentRequest:
    experiment_id: UUID


@dataclass
class LoadConversationRequest:
    experiment_id: UUID


@dataclass
class LoadConversationResponse:
    messages: List[Message]


@dataclass
class CreateMessageRequest:
    experiment_id: UUID
    message_id: UUID
    message_content: str
    role: str


@dataclass
class CreateFunctionCallMessageRequest:
    experiment_id: UUID
    message_id: UUID
    function_call: FunctionCall
    role: str


@dataclass
class CreateFunctionCallMessageResponse:
    saved_message: Message


@dataclass
class CreateMessageResponse:
    saved_message: Message


@dataclass
class SendQueryRequest:
    experiment_id: UUID
    query: str


@dataclass
class SendQueryResponse:
    biobuddy_response: Message


@dataclass
class GetAvailableFunctionCallsResponse:
    function_calls: List[Dict[str, str | Dict[str, dict[str, str] | Any]]]


@dataclass
class EditMessageRequest:
    experiment_id: UUID
    message_id: UUID
    message_content: str


@dataclass
class EditMessageResponse:
    edited_message: Message
