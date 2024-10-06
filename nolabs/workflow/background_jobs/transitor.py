from nolabs.infrastructure.cel import cel
from nolabs.infrastructure.red import redis_client
from workflow.logic.correlation import _all_correlation_ids, _unpack_correlation_id, _all_celery_tasks
from workflow.logic.states import _get_internal_state, ControlStates, _set_internal_state
import celery.states as celery_states


_matching_states = {
    celery_states.FAILURE : (ControlStates.FAILURE),
    celery_states.SUCCESS : (ControlStates.SUCCESS, ControlStates.CANCELLED),
    celery_states.STARTED: (ControlStates.STARTED),
    celery_states.REVOKED: (ControlStates.FAILURE),
    celery_states.RETRY: (ControlStates.FAILURE),
    celery_states.RECEIVED: (ControlStates.PENDING),
    celery_states.PENDING: (ControlStates.PENDING)
}


async def observe_states():
    for correlation_id in await _all_correlation_ids():
        experiment_id, component_id, job_id = _unpack_correlation_id(correlation_id=correlation_id)
        celery_task_ids = await _all_celery_tasks(cid=correlation_id)

        if not job_id:
            # get component task states
            for celery_task_id in celery_task_ids:
                task_result = cel.task_result(task_id=celery_task_id)
                internal_state, internal_state_message = _get_internal_state(id=component_id)

                # states are not matching
                if internal_state not in _matching_states[task_result.state]:
                    internal_state = _matching_states[task_result.state][0]
                    await _set_internal_state(id=component_id, state=internal_state, state_message=task_result.state.message)
        else:


