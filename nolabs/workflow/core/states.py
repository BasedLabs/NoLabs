from enum import Enum

import celery.states as celery_states


class ControlStates(str, Enum):
    SCHEDULED = "SCHEDULED"
    FAILURE = "FAILURE"
    STARTED = "STARTED"
    CANCELLING = "CANCELLING"
    SUCCESS = "SUCCESS"
    CANCELLED = "CANCELLED"
    UNKNOWN = "UNKNOWN"


celery_to_internal_mapping = {
    celery_states.FAILURE: ControlStates.FAILURE,
    celery_states.RETRY: ControlStates.FAILURE,
    celery_states.REVOKED: ControlStates.CANCELLED,
    celery_states.PENDING: ControlStates.STARTED,
    celery_states.RECEIVED: ControlStates.STARTED,
    celery_states.STARTED: ControlStates.STARTED,
    celery_states.SUCCESS: ControlStates.SUCCESS,
    celery_states.REJECTED: ControlStates.FAILURE,
    celery_states.IGNORED: ControlStates.FAILURE,
}

TERMINAL_STATES = [
    ControlStates.FAILURE,
    ControlStates.SUCCESS,
    ControlStates.CANCELLED,
]
ERROR_STATES = [ControlStates.FAILURE]
PROGRESS_STATES = [ControlStates.STARTED, ControlStates.CANCELLING]

state_transitions = {
    ControlStates.UNKNOWN: [
        ControlStates.UNKNOWN,
        ControlStates.SCHEDULED,
        ControlStates.CANCELLING,
        ControlStates.CANCELLED,
        ControlStates.SUCCESS,
    ],
    ControlStates.CANCELLING: [ControlStates.CANCELLED],
    ControlStates.SCHEDULED: [
        ControlStates.SCHEDULED,
        ControlStates.STARTED,
        ControlStates.CANCELLED,
        ControlStates.SUCCESS,
        ControlStates.CANCELLING,
    ],
    ControlStates.STARTED: [
        ControlStates.CANCELLING,
        ControlStates.FAILURE,
        ControlStates.SUCCESS,
        ControlStates.CANCELLED,
    ],
    ControlStates.FAILURE: [
        ControlStates.FAILURE,
        ControlStates.CANCELLING,
        ControlStates.UNKNOWN,
    ],
    ControlStates.CANCELLED: [ControlStates.CANCELLED, ControlStates.UNKNOWN],
    ControlStates.SUCCESS: [
        ControlStates.SUCCESS,
        ControlStates.CANCELLED,
        ControlStates.UNKNOWN,
    ],
}
