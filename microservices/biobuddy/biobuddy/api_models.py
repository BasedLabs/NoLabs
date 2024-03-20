from __future__ import annotations

from pydantic import dataclasses as pcdataclass
from typing import List, Optional, Dict, Any

from openai.types.chat import ChatCompletionMessage

from biobuddy.mixins import BaseModelMixin


@pcdataclass.dataclass
class SendMessageToBioBuddyRequest(BaseModelMixin):
    message_content: str
    # TODO: send only the last N tokens based on model's context window
    previous_messages: List[Dict[str, str]]  # Need to pass them all the time for preserving context
    tools: List[Dict[str, Any]]
    job_id: Optional[str] = None


@pcdataclass.dataclass
class SendMessageToBioBuddyResponse(BaseModelMixin):
    chatgpt_reply: ChatCompletionMessage

@pcdataclass.dataclass
class IsJobRunningResponse:
    is_running: bool
