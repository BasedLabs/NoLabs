import json
import uuid
from enum import Enum
from typing import Optional, Tuple

import celery.states as celery_states

from nolabs.infrastructure.red import redis_client


class ControlStates(str, Enum):
    FAILURE = "FAILURE"
    STARTED = "STARTED"
    CANCELLED = "CANCELLED"
    SUCCESS = "SUCCESS"
    UNKNOWN = "UNKNOWN"


UNREADY_STATES = [ControlStates.STARTED, ControlStates.UNKNOWN]
TERMINAL_STATES = [ControlStates.FAILURE, ControlStates.SUCCESS, ControlStates.CANCELLED]


celery_to_internal_mapping = {
    celery_states.FAILURE: ControlStates.FAILURE,
    celery_states.RETRY:  ControlStates.FAILURE,
    celery_states.REVOKED:  ControlStates.CANCELLED,
    celery_states.PENDING: ControlStates.STARTED,
    celery_states.RECEIVED: ControlStates.STARTED,
    celery_states.STARTED: ControlStates.STARTED,
    celery_states.SUCCESS: ControlStates.SUCCESS,
    celery_states.REJECTED: ControlStates.CANCELLED,
    celery_states.IGNORED: ControlStates.UNKNOWN
}


async def _set_internal_state(id: uuid.UUID | str, state: ControlStates, state_message: Optional[str] = None):
    key = f"state:{id}"
    j = json.dumps({'state': state, 'state_message': state_message})
    await redis_client.set(name=key, value=j)


async def _get_internal_state(id: uuid.UUID | str) -> Tuple[ControlStates, Optional[str]]:
    key = f"state:{id}"
    obj = await redis_client.get(key)

    if obj is None:
        return ControlStates.UNKNOWN, None

    obj = json.loads(obj)

    return ControlStates(obj['state']), obj['state_message']

async def _ready(id: uuid.UUID) -> bool:
    state, _ = await _get_internal_state(id)
    return state not in UNREADY_STATES