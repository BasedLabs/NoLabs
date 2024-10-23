import uuid
from typing import List, Iterable

from nolabs.domain.models.common import Job
from nolabs.infrastructure.redis_client_factory import get_redis_pipe
from workflow.core import Tasks
from nolabs.workflow.core.job_execution_nodes import JobExecutionNode
from nolabs.workflow.core.node import CeleryExecutionNode, ExecutionNode
from nolabs.workflow.core.socketio_events_emitter import emit_start_component_event
from nolabs.workflow.core.states import ControlStates
from nolabs.workflow.core.states import TERMINAL_STATES


class ComponentMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self,
                 experiment_id: uuid.UUID,
                 component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:main_task")
        self.experiment_id = experiment_id
        self.component_id = component_id

    async def start(self):
        pipe = get_redis_pipe()
        celery_task_id = await self._prepare_for_start(pipe=pipe)
        self.celery.send_task(name=Tasks.component_main_task,
                              kwargs={"experiment_id": self.experiment_id,
                                      "component_id": self.component_id},
                              task_id=celery_task_id,
                              retry=False
                              )
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        await self.set_output(output={}, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        await pipe.execute()

    async def schedule(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        await self.set_state(state=ControlStates.SCHEDULED)


class ComponentCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:complete_task")
        self.experiment_id = experiment_id
        self.component_id = component_id

    async def start(self):
        pipe = get_redis_pipe()
        celery_task_id = await self._prepare_for_start(pipe=pipe)
        self.celery.send_task(name=Tasks.complete_component_task,
                              kwargs={"experiment_id": self.experiment_id,
                                      "component_id": self.component_id},
                              task_id=celery_task_id,
                              retry=False
                              )
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        await self.set_output(output={}, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        await pipe.execute()

    async def schedule(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_ids: List[uuid.UUID]):
        await self.set_state(state=ControlStates.SCHEDULED)


class ComponentExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}")
        self.component_id = component_id
        self.experiment_id = experiment_id
        self.main_task = ComponentMainTaskExecutionNode(experiment_id, component_id)
        self.complete_task = ComponentCompleteTaskExecutionNode(experiment_id, component_id)

    async def can_schedule(self, previous_component_ids: Iterable[uuid.UUID]) -> bool:
        if not await super().can_schedule():
            return False

        for component_id in previous_component_ids:
            component_node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component_id)
            if await component_node.get_state() != ControlStates.SUCCESS:
                return False

        return True

    async def start(self, **kwargs):
        await self.set_state(state=ControlStates.STARTED)
        emit_start_component_event(experiment_id=self.experiment_id,
                                   component_id=self.component_id)

    async def sync_started(self):
        state = await self.get_state()
        if state != ControlStates.STARTED:
            return

        await self.sync_main_task()

        if await self.get_state() in TERMINAL_STATES:
            return

        any_job_succeeded = False
        if await self.main_task.get_state() == ControlStates.SUCCESS:
            any_job_succeeded = await self.sync_jobs()
            if await self.get_state() in TERMINAL_STATES:
                return

        if await self.get_state() in TERMINAL_STATES:
            return

        if await self.main_task.get_state() == ControlStates.SUCCESS and any_job_succeeded:
            await self.sync_complete_task()

        if await self.get_state() in TERMINAL_STATES:
            return

    async def sync_main_task(self):
        if await self.main_task.can_schedule():
            await self.main_task.schedule(
                experiment_id=self.experiment_id,
                component_id=self.component_id
            )

        if await self.main_task.can_start():
            await self.main_task.start()
            emit_start_component_event(experiment_id=self.experiment_id, component_id=self.component_id)

        await self.main_task.sync_started()

        state = await self.main_task.get_state()
        if state in TERMINAL_STATES and state != ControlStates.SUCCESS:
            await self.set_state(await self.main_task.get_state())
            await self.set_message(await self.main_task.get_message())

    async def sync_jobs(self, batch_size=4) -> bool:
        """
        :returns: any jobs succeeded
        """
        job_ids = [j.id for j in Job.objects(component=self.component_id).only('id')]
        active_jobs = 0

        # Track job nodes
        for job_id in job_ids:
            job_node = JobExecutionNode(self.experiment_id, self.component_id, job_id)
            await job_node.sync_started()

            job_state = await job_node.get_state()
            if job_state in TERMINAL_STATES:
                continue  # Skip jobs in terminal state

            if active_jobs >= batch_size:
                return False  # Stop tracking more jobs

            active_jobs += 1

            if await job_node.can_schedule():
                await job_node.schedule(
                    experiment_id=self.experiment_id,
                    component_id=self.component_id,
                    job_id=job_id
                )
            elif await job_node.can_start():
                await job_node.start()

        # If no jobs are running and main task is completed, update component state
        if active_jobs == 0 and await self.main_task.get_state() in TERMINAL_STATES:
            return await self.update_job_states(job_ids)

        return False

    async def update_job_states(self, job_ids: List[uuid.UUID]) -> bool:
        """
        Update the component state based on the terminal states of all jobs.
        :returns: any jobs succeeded
        """

        if not job_ids:
            return True

        all_jobs_terminal = True
        any_success = False

        for job_id in job_ids:
            job_node = JobExecutionNode(self.experiment_id, self.component_id, job_id)
            job_state = await job_node.get_state()

            if job_state in TERMINAL_STATES:
                if job_state == ControlStates.SUCCESS:
                    any_success = True
            else:
                all_jobs_terminal = False

        # If all jobs are in terminal state, update component state
        if all_jobs_terminal:
            if not any_success:
                pipe = get_redis_pipe()
                await self.set_state(ControlStates.FAILURE, pipe=pipe)
                await self.set_message(message="All jobs failed", pipe=pipe)
                await pipe.execute()
                return any_success

        return any_success

    async def sync_complete_task(self):
        """Synchronize and start the complete task if all jobs are finished."""
        if await self.complete_task.get_state() in TERMINAL_STATES:
            return

        await self.complete_task.sync_started()

        job_ids = [j.id for j in Job.objects(component=self.component_id).only('id')]

        # If all jobs are in terminal states and any is successful, run the complete task
        if await self.complete_task.can_schedule():
            await self.complete_task.schedule(
                experiment_id=self.experiment_id,
                component_id=self.component_id,
                job_ids=job_ids
            )

        if await self.complete_task.can_start():
            await self.complete_task.start()

        # Update the component state if the complete task finished
        if await self.complete_task.get_state() in TERMINAL_STATES:
            await self.set_state(await self.complete_task.get_state())
            await self.set_message(await self.complete_task.get_message())
