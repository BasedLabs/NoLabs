import uuid
from typing import Iterable, List

from redis.client import Pipeline

from nolabs.domain.models.common import Job
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.redis_client_factory import get_redis_pipe
from nolabs.infrastructure.settings import settings
from nolabs.workflow.core import Tasks
from nolabs.workflow.core.job_execution_nodes import JobExecutionNode
from nolabs.workflow.core.node import CeleryExecutionNode, ExecutionNode
from nolabs.workflow.core.socketio_events_emitter import (
    emit_component_jobs_event,
    emit_finish_component_event,
    emit_start_component_event,
)
from nolabs.workflow.core.states import (
    ERROR_STATES,
    PROGRESS_STATES,
    TERMINAL_STATES,
    ControlStates,
)


class ComponentMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:main_task")
        self.experiment_id = experiment_id
        self.component_id = component_id

    async def start(self):
        logger.debug(
            "Starting component main task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            }
        )
        pipe = get_redis_pipe()
        queue = settings.workflow_queue
        celery_task_id = await self._prepare_for_start(queue=queue, pipe=pipe)
        self.celery.send_task(
            name=Tasks.component_main_task,
            kwargs={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
            queue=queue,
            task_id=celery_task_id,
            retry=False,
        )
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        pipe.execute()

    async def schedule(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        logger.debug(
            "Scheduling component main task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
        )
        await self.set_state(state=ControlStates.SCHEDULED)

    async def on_finished(self):
        await emit_component_jobs_event(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_ids=[j.id for j in Job.objects(component=self.component_id).only("id")],
        )


class ComponentCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(
            id=f"execution_node:{experiment_id}:{component_id}:complete_task"
        )
        self.experiment_id = experiment_id
        self.component_id = component_id

    async def start(self):
        logger.debug(
            "Starting component complete task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
        )
        pipe = get_redis_pipe()
        queue = settings.workflow_queue
        celery_task_id = await self._prepare_for_start(queue=queue, pipe=pipe)
        self.celery.send_task(
            name=Tasks.complete_component_task,
            kwargs={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
            queue=queue,
            task_id=celery_task_id,
            retry=False,
        )
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        pipe.execute()

    async def schedule(
        self,
        experiment_id: uuid.UUID,
        component_id: uuid.UUID,
        job_ids: List[uuid.UUID],
    ):
        logger.debug(
            "Scheduling component complete task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
        )
        await self.set_state(state=ControlStates.SCHEDULED)


class ComponentExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}")
        self.component_id = component_id
        self.experiment_id = experiment_id
        self.main_task = ComponentMainTaskExecutionNode(experiment_id, component_id)
        self.complete_task = ComponentCompleteTaskExecutionNode(
            experiment_id, component_id
        )

    async def can_start(self, previous_component_ids: Iterable[uuid.UUID]) -> bool:
        if not await super().can_start():
            return False

        for component_id in previous_component_ids:
            component_node = ComponentExecutionNode(
                experiment_id=self.experiment_id, component_id=component_id
            )
            if await component_node.get_state() != ControlStates.SUCCESS:
                return False

        return True

    async def can_schedule(self) -> bool:
        state = await self.get_state()
        return state == ControlStates.UNKNOWN or state in TERMINAL_STATES

    async def schedule(self):
        logger.debug(
            "Scheduling component node",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
        )
        pipe = get_redis_pipe()
        await self.reset(pipe=pipe)
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        pipe.execute()

    async def reset(self, pipe: Pipeline):
        state = await self.get_state()

        if state not in TERMINAL_STATES:
            return

        await super().reset(pipe=pipe)
        main_task = ComponentMainTaskExecutionNode(
            experiment_id=self.experiment_id, component_id=self.component_id
        )
        await main_task.reset(pipe=pipe)
        job_ids = [j.id for j in Job.objects(component=self.component_id).only("id")]
        for job_id in job_ids:
            job_node = JobExecutionNode(self.experiment_id, self.component_id, job_id)
            job_state = await job_node.get_state()
            if job_state in ERROR_STATES and job_state != ControlStates.SUCCESS:
                await job_node.reset(pipe=pipe)
        complete_task = ComponentCompleteTaskExecutionNode(
            experiment_id=self.experiment_id, component_id=self.component_id
        )
        await complete_task.reset(pipe=pipe)

    async def start(self, **kwargs):
        logger.debug(
            "Starting component node",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
        )
        await self.set_state(state=ControlStates.STARTED)

    async def sync_started(self):
        state = await self.get_state()
        logger.debug(
            "Syncing component node",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "state": str(state),
            },
        )
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

        if (
            await self.main_task.get_state() == ControlStates.SUCCESS
            and any_job_succeeded
        ):
            await self.sync_complete_task()

        if await self.get_state() in TERMINAL_STATES:
            return

    async def sync_main_task(self):
        logger.debug(
            "Syncing component main task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
        )
        if await self.main_task.can_schedule():
            await self.main_task.schedule(
                experiment_id=self.experiment_id, component_id=self.component_id
            )

        if await self.main_task.can_start():
            await self.main_task.start()

        await self.main_task.sync_started()

        state = await self.main_task.get_state()
        if state in TERMINAL_STATES and state != ControlStates.SUCCESS:
            await self.set_state(await self.main_task.get_state())
            await self.set_message(await self.main_task.get_message())

    async def sync_jobs(self, batch_size=4) -> bool:
        """
        :returns: any jobs succeeded
        """
        job_ids = [j.id for j in Job.objects(component=self.component_id).only("id")]
        logger.debug(
            "Syncing component jobs",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "jobs_count": len(job_ids),
            },
        )
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
                await job_node.schedule()
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
        logger.debug(
            "Updating component jobs states",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "jobs_count": len(job_ids),
            },
        )

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
                pipe.execute()
                return any_success

        return any_success

    async def sync_complete_task(self):
        """Synchronize and start the complete task if all jobs are finished."""
        logger.debug(
            "Syncing component complete task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
        )
        if await self.complete_task.get_state() in TERMINAL_STATES:
            await self.set_state(await self.complete_task.get_state())
            await self.set_message(await self.complete_task.get_message())
            return

        await self.complete_task.sync_started()

        job_ids = [j.id for j in Job.objects(component=self.component_id).only("id")]

        # If all jobs are in terminal states and any is successful, run the complete task
        if await self.complete_task.can_schedule():
            await self.complete_task.schedule(
                experiment_id=self.experiment_id,
                component_id=self.component_id,
                job_ids=job_ids,
            )

        if await self.complete_task.can_start():
            await self.complete_task.start()

    async def on_started(self):
        await emit_start_component_event(
            experiment_id=self.experiment_id, component_id=self.component_id
        )

    async def on_finished(self):
        await emit_finish_component_event(
            experiment_id=self.experiment_id, component_id=self.component_id
        )

    async def sync_cancelling(self):
        logger.debug(
            "Cancelling component",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
            },
        )
        if await self.get_state() != ControlStates.CANCELLING:
            return

        if await self.main_task.can_cancel():
            await self.main_task.cancel()

        await self.main_task.sync_cancelling()

        jobs_in_progress = True
        job_ids = [j.id for j in Job.objects(component=self.component_id).only("id")]
        for job_id in job_ids:
            job_node = JobExecutionNode(self.experiment_id, self.component_id, job_id)
            if await job_node.can_cancel():
                await job_node.cancel()
            jobs_in_progress = await job_node.get_state() not in PROGRESS_STATES

        if not jobs_in_progress:
            return

        if await self.complete_task.can_cancel():
            await self.complete_task.cancel()

        await self.complete_task.sync_cancelling()

        if (
            await self.main_task.get_state() not in PROGRESS_STATES
            and await self.complete_task.get_state() not in PROGRESS_STATES
        ):
            await self.set_cancelled()
            return ControlStates.CANCELLED

        return await self.get_state()
