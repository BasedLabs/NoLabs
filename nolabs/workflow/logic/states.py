import json
import uuid
from enum import Enum
from typing import Optional, Tuple

from celery.result import AsyncResult
from celery.states import FAILURE, PENDING, STARTED, RETRY, RECEIVED, SUCCESS

from nolabs.infrastructure.red import redis_client


class ControlStates(str, Enum):
    FAILURE = "FAILURE"
    CANCELLED = "CANCELLED"
    PENDING = "PENDING"
    STARTED = "STARTED"
    SUCCESS = "SUCCESS"
    UNKNOWN = "UNKNOWN"


UNREADY_STATES = [ControlStates.PENDING, ControlStates.STARTED]
TERMINAL_STATES = [ControlStates.FAILURE, ControlStates.CANCELLED, ControlStates.SUCCESS]


def _celery_to_internal(result: AsyncResult) -> ControlStates:
    if result.state == FAILURE:
        return ControlStates.FAILURE
    if result.state in [PENDING, RETRY, RECEIVED]:
        return ControlStates.PENDING
    if result.state == STARTED:
        return ControlStates.STARTED
    if result.state == SUCCESS:
        return ControlStates.SUCCESS

    return ControlStates.UNKNOWN


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