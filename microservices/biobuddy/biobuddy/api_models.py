from __future__ import annotations

import dataclasses
from typing import List, Optional

from biobuddy.mixins import BaseModelMixin


@dataclasses.dataclass
class SendMessageToBioBuddyRequest(BaseModelMixin):
    message_content: str
    # TODO: send only the last N tokens based on model's context window
    previous_messages: List[str]  # Need to pass them all the time for preserving context
    job_id: Optional[str] = None


@dataclasses.dataclass
class SendMessageToBioBuddyResponse(BaseModelMixin):
    chatgpt_reply: str

@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool
