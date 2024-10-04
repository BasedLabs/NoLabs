import uuid
from enum import Enum
from typing import Optional, Tuple

from nolabs.infrastructure.red import redis_client


class ControlStates(str, Enum):
    FAILURE = "FAILURE"
    PENDING = "PENDING"
    STARTED = "STARTED"
    SUCCESS = "SUCCESS"


def _set_state(id: uuid.UUID, state: ControlStates, state_message: Optional[str] = None):
    key = f"{id}:state"
    redis_client.set(name=key, value=state)

    if state_message:
        key = f"{key}:message"
        redis_client.set(name=key, value=state_message)


def _get_state(id: uuid.UUID) -> Tuple[str, Optional[str]]:
    key = f"{id}:state"
    state = redis_client.get(key)
    state_message = None

    if state_message:
        key = f"{key}:message"
        state_message = redis_client.get(key)

    return state, state_message
