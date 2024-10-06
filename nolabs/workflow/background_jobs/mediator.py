import uuid
from typing import List

from celery.states import READY_STATES

from domain.models.common import Job
from nolabs.infrastructure.cel import cel
from workflow.logic.correlation import _all_components_cids, _all_celery_tasks, unpack_cid, _make_job_cid, \
    _clear_correlation, _make_component_cid
from workflow.logic.states import _set_internal_state, celery_to_internal_mapping, ControlStates, TERMINAL_STATES, \
    _get_internal_state


@cel.app.task(name="get_executing_components")
async def get_executing_components():
    chunk_size = 10
    cids = await _all_components_cids()
    return [cids[i:i + chunk_size] for i in range(0, len(cids), chunk_size)]


@cel.app.task(name="state_transition_doer")
async def state_transition_doer(component_cids: List[str]):
    for cid in component_cids:
        component_id = unpack_cid(cid=cid)
        job_states = await _transit_job_states(component_id=component_id)
        await _transit_component_state(component_id=component_id, job_states=job_states)

async def _transit_component_state(component_id: uuid.UUID, job_states: List[ControlStates]):
    cid = _make_component_cid(id=component_id)
    celery_task_ids = await _all_celery_tasks(cid=cid)
    if celery_task_ids:
        celery_task_id = celery_task_ids[-1]
        task_result = cel.task_result(task_id=celery_task_id)
        internal_state = celery_to_internal_mapping[task_result.state]

        jobs_are_terminal = not job_states or all(job_state in TERMINAL_STATES for job_state in job_states)
        component_is_terminal = internal_state in TERMINAL_STATES

        if jobs_are_terminal and not component_is_terminal:
            # Jobs are completed, but component is pending


        if jobs_are_terminal and component_is_terminal:
            await _set_internal_state(id=component_id, state=internal_state, state_message=task_result.result)
            await _clear_correlation(cid=cid)

async def _sync_jobs_states(component_id: uuid.UUID) -> List[ControlStates]:
    """
    Synchronize internal jobs states with actual worker state in case of worker process failure
    """
    jobs = Job.objects(component=component_id).get('id')

    job_states: List[ControlStates] = []

    for job in jobs:
        cid = _make_job_cid(id=job.id)
        internal_state, message = _get_internal_state(id=job.id)

        # If already in terminal state (calling code completed the job)
        if internal_state in TERMINAL_STATES:
            await _clear_correlation(cid=cid)
            job_states.append(internal_state)
            continue

        celery_task_ids = await _all_celery_tasks(cid=cid)
        if celery_task_ids:
            celery_task_id = celery_task_ids[-1]
            task_result = cel.task_result(task_id=celery_task_id)
            internal_state = celery_to_internal_mapping[task_result.state]
            if internal_state == ControlStates.FAILURE:
                await _set_internal_state(id=job.id, state=internal_state, state_message=task_result.result)
                await _clear_correlation(cid=cid)
                job_states.append(internal_state)

    return job_states
