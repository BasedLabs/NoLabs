from pydantic import BaseModel
from typing import Optional, List
from enum import Enum

class MessageType(str, Enum):
    human = 'human'
    ai = 'ai'
    tool = 'tool'

class ChatRequest(BaseModel):
    message: str

class DocumentContext(BaseModel):
    page_content: str
    title: str
    id: str

class ChatResponse(BaseModel):
    id: str
    context: List[DocumentContext]
    message_type: MessageType
    content: Optional[str] = None

class SetOpenAiApiKey(BaseModel):
    api_key: str
