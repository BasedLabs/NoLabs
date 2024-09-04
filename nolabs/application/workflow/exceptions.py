from enum import Enum
from typing import List

from pydantic.dataclasses import dataclass


@dataclass
class ErrorCode:
    """
    Will be used to represent and describe the error in application
    """
    code: int
    description: str

class ErrorCodes(Enum):
    invalid_task_state = ErrorCode(code=1, description='Invalid task state')


if len([e.value.code for e in ErrorCodes]) != len(set([e.value.code for e in ErrorCodes])):
    raise ValueError("Invalid ErrorCode initialization")


class WorkflowException(Exception):
    def __init__(self, error_code: ErrorCodes, messages: str | List[str] | None = None):
        self.error_code = error_code.value.code

        if messages:
            if isinstance(messages, str):
                self.messages = [messages]
            else:
                self.messages = messages
        else:
            self.messages = [
                error_code.value.description
            ]