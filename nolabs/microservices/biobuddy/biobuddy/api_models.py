from __future__ import annotations

from typing import Any, Dict, List, Optional

from biobuddy.mixins import BaseModelMixin
from pydantic import dataclasses as pcdataclass


@pcdataclass.dataclass
class Component:
    name: str
    description: str | None = None
    inputs = List[str]
    outputs = List[str]


@pcdataclass.dataclass
class Connection:
    source_path: List[str]
    target_path: List[str]
    source_component_id: str


@pcdataclass.dataclass
class WorkflowComponent:
    id: str
    name: str
    description: str
    connections: List[Connection]


@pcdataclass.dataclass
class SendMessageToBioBuddyRequest(BaseModelMixin):
    experiment_id: str
    message_content: str
    # TODO: send only the last N tokens based on model's context window
    previous_messages: List[Dict[str, str]]
    available_components: List[Component]
    tools: List[Dict[str, Any]]
    current_workflow: List[WorkflowComponent] = None
    job_id: Optional[str] = None


@pcdataclass.dataclass
class SendActionCallRequest(BaseModelMixin):
    action_text: str
    tools: List[Dict[str, Any]]
    job_id: Optional[str] = None


@pcdataclass.dataclass
class SendMessageToBioBuddyResponse(BaseModelMixin):
    reply_type: str
    content: str

    def to_dict(self):
        return {
            "reply_type": self.reply_type,
            "content": self.content,
        }


@pcdataclass.dataclass
class IsJobRunningResponse:
    is_running: bool
