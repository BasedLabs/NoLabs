import uuid
from typing import List

from domain.models.common import Job
from infrastructure.log import get_worker_logger
from infrastructure.redis_client_factory import use_redis_pipe
from nolabs.workflow.core.celery_tasks import _component_main_task, _complete_component_task
from nolabs.workflow.core.job_execution_nodes import JobExecutionNode
from nolabs.workflow.core.node import CeleryExecutionNode, ExecutionNode
from nolabs.workflow.core.socketio_events_emitter import emit_start_component_event
from nolabs.workflow.core.states import ControlStates
from nolabs.workflow.core.states import TERMINAL_STATES

logger = get_worker_logger()


class ComponentMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self,
                 experiment_id: uuid.UUID,
                 component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:main_task")
        self.experiment_id = experiment_id
        self.component_id = component_id

    async def execute(self):
        async with use_redis_pipe():
            celery_task_id = await self._prepare_for_execute()
            _component_main_task.apply_async(
                kwargs={"experiment_id": self.experiment_id,
                        "component_id": self.component_id},
                task_id=celery_task_id,
                retry=False)
            await self.set_state(ControlStates.STARTED)
            await self.set_output(output={})
            await self.set_message(message="")

    async def schedule(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        async with use_redis_pipe():
            await self.set_state(state=ControlStates.SCHEDULED)


class ComponentCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:complete_task")
        self.experiment_id = experiment_id
        self.component_id = component_id

    async def can_schedule(self, job_ids: List[uuid.UUID]):
        states = []
        for job_id in job_ids:
            job_node = JobExecutionNode(experiment_id=self.experiment_id, component_id=self.component_id, job_id=job_id)
            job_state = await job_node.get_state()
            if job_state not in TERMINAL_STATES:
                return False
            states.append(job_state)
        return await super().can_schedule() and any(s == ControlStates.SUCCESS for s in states)

    async def execute(self):
        async with use_redis_pipe():
            celery_task_id = await self._prepare_for_execute()
            _complete_component_task.apply_async(
                kwargs={"experiment_id": self.experiment_id,
                        "component_id": self.component_id},
                task_id=celery_task_id,
                retry=False)
            await self.set_state(ControlStates.STARTED)
            await self.set_output(output={})
            await self.set_message(message="")

    async def schedule(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_ids: List[uuid.UUID]):
        async with use_redis_pipe():
            await self.set_state(state=ControlStates.SCHEDULED)


class ComponentExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}")
        self.component_id = component_id
        self.experiment_id = experiment_id
        self.main_task = ComponentMainTaskExecutionNode(experiment_id, component_id)
        self.complete_task = ComponentCompleteTaskExecutionNode(experiment_id, component_id)

    async def can_execute(self, previous_component_ids: List[uuid.UUID]):
        if not super().can_execute():
            return False

        for component_id in previous_component_ids:
            component_node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component_id)
            if await component_node.get_state() != ControlStates.SUCCESS:
                return False

        return True

    async def execute(self, **kwargs):
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

        if await self.main_task.get_state() == ControlStates.SUCCESS and any_job_succeeded:
            await self.sync_complete_task()

    async def sync_main_task(self):
        await self.main_task.sync_started()

        if await self.main_task.can_schedule():
            await self.main_task.schedule(
                experiment_id=self.experiment_id,
                component_id=self.component_id
            )

        if await self.main_task.can_execute():
            await self.main_task.execute()
            emit_start_component_event(experiment_id=self.experiment_id, component_id=self.component_id)

        if await self.main_task.get_state() in TERMINAL_STATES:
            await self.set_state(await self.main_task.get_state())
            await self.set_message(await self.main_task.get_message())

    async def sync_jobs(self, batch_size=4) -> bool:
        """
        :returns: any jobs succeeded
        """
        input_data = (await self.get_input()) or dict()
        job_ids = input_data.get('job_ids', [])
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
                    experiment_id=input_data["experiment_id"],
                    component_id=input_data["component_id"],
                    job_id=job_id
                )
            elif await job_node.can_execute():
                await job_node.execute()

        # If no jobs are running and main task is completed, update component state
        if active_jobs == 0 and await self.main_task.get_state() in TERMINAL_STATES:
            return await self.update_job_states(job_ids)

        return False

    async def update_job_states(self, job_ids: List[uuid.UUID]) -> bool:
        """
        Update the component state based on the terminal states of all jobs.
        :returns: any jobs succeeded
        """
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
            if any_success:
                await self.set_state(ControlStates.SUCCESS)
            else:
                await self.set_state(ControlStates.FAILURE)

        return any_success

    async def sync_complete_task(self):
        """Synchronize and execute the complete task if all jobs are finished."""
        if await self.complete_task.get_state() in TERMINAL_STATES:
            return

        job_ids = [j.id for j in Job.objects(component=self.component_id).only('id')]

        # If all jobs are in terminal states and any is successful, run the complete task
        if await self.complete_task.can_schedule(job_ids):
            await self.complete_task.schedule(
                experiment_id=self.experiment_id,
                component_id=self.component_id,
                job_ids=job_ids
            )

        if await self.complete_task.can_execute():
            await self.complete_task.execute()

        # Update the component state if the complete task finished
        if await self.complete_task.get_state() in TERMINAL_STATES:
            await self.set_state(await self.complete_task.get_state())
            await self.set_message(await self.complete_task.get_message())
