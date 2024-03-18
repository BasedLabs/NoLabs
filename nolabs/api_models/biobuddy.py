from pydantic import dataclasses as pcdataclass, BaseModel
from typing import List, Any, Union, Optional, Dict


@pcdataclass.dataclass
class CheckBioBuddyEnabledResponse:
    enabled: bool


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


class FileData(BaseModel):
    content: Optional[str] = None  # Assuming file content is a string; adjust as needed
    metadata: Dict[str, Any]  # Metadata related to the file


# Subclass for RcsbPdbData
@pcdataclass.dataclass
class RcsbPdbMetaData:
    link: str

class RcsbPdbData(FileData):
    metadata: RcsbPdbMetaData  # Example of a RcsbPdbMetaData

# Subclass for ChemBLData
@pcdataclass.dataclass
class ChemBLMetaData:
    chembl_id: str
    link: str
    pref_name: str


class ChemBLData(FileData):
    metadata: ChemBLMetaData  # Example of a specific attribute for ChemBLData

# Base class for return data, which now contains a list of FileData objects
class FunctionCallReturnData(BaseModel):
    files: List[Union[
        RcsbPdbData,
        ChemBLData
    ]]  # List of files with their content and metadata


@pcdataclass.dataclass
class FunctionCall:
    function_name: str
    parameters: List[FunctionParam]
    data: FunctionCallReturnData = None

@pcdataclass.dataclass
class Message:
    role: str
    message: Union[RegularMessage, FunctionCall]
    type: str


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
    biobuddy_response: Message
