import uuid
from typing import Any, Dict, Optional

from celery.result import AsyncResult, allow_join_result
from redis.client import Pipeline

from nolabs.infrastructure.log import logger
from nolabs.infrastructure.redis_client_factory import get_redis_pipe, rd
from nolabs.infrastructure.settings import settings
from nolabs.workflow.core import Tasks
from nolabs.workflow.core.node import CeleryExecutionNode, ExecutionNode
from nolabs.workflow.core.socketio_events_emitter import (
    emit_finish_job_event,
    emit_start_job_event,
)
from nolabs.workflow.core.states import (
    ERROR_STATES,
    PROGRESS_STATES,
    TERMINAL_STATES,
    ControlStates,
)


class JobMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(
        self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID
    ):
        super().__init__(
            id=f"execution_node:{experiment_id}:{component_id}:{job_id}:main_task"
        )
        self.experiment_id = experiment_id
        self.component_id = component_id
        self.job_id = job_id

    async def start(self):
        logger.debug(
            "Starting job main task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )
        pipe = get_redis_pipe()
        queue = settings.workflow_queue
        celery_task_id = await self._prepare_for_start(queue=queue, pipe=pipe)
        self.celery.send_task(
            name=Tasks.job_main_task,
            kwargs={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
            task_id=celery_task_id,
            retry=False,
        )
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        pipe.execute()

    async def schedule(self):
        logger.debug(
            "Scheduling job main task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )
        pipe = get_redis_pipe()
        await self.set_input(
            execution_input={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
            pipe=pipe,
        )
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        pipe.execute()


class JobLongRunningTaskExecutionNode(CeleryExecutionNode):
    def __init__(
        self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID
    ):
        super().__init__(
            id=f"execution_node:{experiment_id}:{component_id}:{job_id}:long_running_task"
        )
        self.experiment_id = experiment_id
        self.component_id = component_id
        self.job_id = job_id

    async def start(self):
        logger.debug(
            "Starting job long running task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )

        input_data = await self.get_input()
        pipe = get_redis_pipe()
        celery_task_name = input_data["celery_task_name"]
        celery_queue = input_data["celery_queue"]
        kwargs = input_data["kwargs"]
        task_id = await self._prepare_for_start(queue=celery_queue, pipe=pipe)
        self.celery.send_task(
            name=celery_task_name,
            task_id=task_id,
            queue=celery_queue,
            retry=False,
            args=[kwargs[list(kwargs.keys())[0]]],
        )
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        pipe.execute()

    async def schedule(
        self,
        celery_task_name: str,
        arguments: Optional[Dict[str, Any]] = None,
        celery_queue: Optional[str] = None,
    ):
        logger.debug(
            "Scheduling job long running task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )
        data = {
            "celery_task_name": celery_task_name,
            "celery_queue": celery_queue,
            "kwargs": arguments or {},
        }

        pipe = get_redis_pipe()
        await self.set_input(execution_input=data, pipe=pipe)
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        pipe.execute()

    async def get_output(self) -> Dict[str, Any]:
        celery_task_id = self._get_task_id()
        if not celery_task_id:
            return {}
        async_result = AsyncResult(id=celery_task_id, app=self.celery)
        if not async_result.ready():
            return {}
        with allow_join_result():
            return async_result.get()


class JobCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(
        self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID
    ):
        super().__init__(
            id=f"execution_node:{experiment_id}:{component_id}:{job_id}:complete_task"
        )
        self.job_id = job_id
        self.component_id = component_id
        self.experiment_id = experiment_id

    async def start(self):
        logger.debug(
            "Starting job complete task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )
        input_data = await self.get_input()
        pipe = get_redis_pipe()
        queue = settings.workflow_queue
        celery_task_id = await self._prepare_for_start(queue=queue, pipe=pipe)
        self.celery.send_task(
            name=Tasks.complete_job_task,
            kwargs={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
                "long_running_output": input_data["long_running_output"],
            },
            task_id=celery_task_id,
            retry=False,
        )
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        pipe.execute()

    async def schedule(self, long_running_output: Optional[Dict[str, Any]] = None):
        logger.debug(
            "Scheduling job complete task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )

        pipe = get_redis_pipe()
        await self.set_input(
            execution_input={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
                "long_running_output": long_running_output or {},
            },
            pipe=pipe,
        )
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        pipe.execute()


class JobExecutionNode(ExecutionNode):
    def __init__(
        self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID
    ):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}")
        self.experiment_id = experiment_id
        self.component_id = component_id
        self.job_id = job_id

        self.main_task = JobMainTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id,
        )

        self.long_running_task = JobLongRunningTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id,
        )

        self.complete_task = JobCompleteTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id,
        )

    async def sync_started(self):
        logger.debug(
            "Syncing job node",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )

        if await self.get_state() != ControlStates.STARTED:
            return

        await self.sync_main_task()
        if await self.get_state() in TERMINAL_STATES:
            return

        await self.sync_long_running_task()
        if await self.get_state() in TERMINAL_STATES:
            return

        await self.sync_complete_task()
        if await self.get_state() in TERMINAL_STATES:
            return

        await self.set_final_state()

    async def sync_main_task(self):
        logger.debug(
            "Syncing job main task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )

        if await self.main_task.can_schedule():
            await self.main_task.schedule()

        if await self.main_task.can_start():
            await self.main_task.start()

        await self.main_task.sync_started()

    async def sync_long_running_task(self):
        logger.debug(
            "Syncing job long running task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )
        main_task_state = await self.main_task.get_state()

        # will not run
        if (
            main_task_state in TERMINAL_STATES
            and await self.long_running_task.get_state() == ControlStates.UNKNOWN
        ):
            await self.long_running_task.set_state(state=ControlStates.SUCCESS)
            return

        if await self.long_running_task.can_start():
            await self.long_running_task.start()
            return

        await self.long_running_task.sync_started()

    async def sync_complete_task(self):
        logger.debug(
            "Syncing job complete task",
            extra={
                "experiment_id": self.experiment_id,
                "component_id": self.component_id,
                "job_id": self.job_id,
            },
        )
        long_running_task_state = await self.long_running_task.get_state()
        main_task_state = await self.main_task.get_state()

        if main_task_state == ControlStates.SUCCESS and long_running_task_state in [
            ControlStates.SUCCESS,
            ControlStates.UNKNOWN,
        ]:
            if await self.complete_task.can_schedule():
                await self.complete_task.schedule(
                    long_running_output=await self.long_running_task.get_output()
                )
                return

        if await self.complete_task.can_start():
            await self.complete_task.start()
            return

        await self.complete_task.sync_started()

    async def set_final_state(self):
        success = True

        for task in [self.main_task, self.long_running_task, self.complete_task]:
            state = await task.get_state()

            if state not in TERMINAL_STATES:
                return

            if state in ERROR_STATES:
                success = False
                await self.set_message(await task.get_message())
                await self.set_state(state)
                break

        if success:
            await self.set_state(ControlStates.SUCCESS)
            await self.set_message("")

    async def start(self):
        await self.set_state(state=ControlStates.STARTED)

    async def schedule(self):
        pipe = get_redis_pipe()
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        await self.reset(pipe=pipe)
        pipe.execute()

    async def reset(self, pipe: Pipeline):
        state = await self.get_state()

        if state not in TERMINAL_STATES:
            return

        await super().reset(pipe=pipe)
        main_task = JobMainTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id,
        )
        await main_task.reset(pipe=pipe)
        long_running_task = JobLongRunningTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id,
        )
        await long_running_task.reset(pipe=pipe)
        complete_task = JobCompleteTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id,
        )
        await complete_task.reset(pipe=pipe)

    async def on_started(self):
        await emit_start_job_event(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id,
        )

    async def on_finished(self):
        await emit_finish_job_event(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id,
        )

    async def sync_cancelling(self):
        if await self.get_state() != ControlStates.CANCELLING:
            return

        if await self.main_task.can_cancel():
            await self.main_task.cancel()

        await self.main_task.sync_cancelling()

        if await self.long_running_task.can_cancel():
            await self.long_running_task.cancel()

        await self.long_running_task.sync_cancelling()

        if await self.complete_task.can_cancel():
            await self.complete_task.cancel()

        await self.complete_task.sync_cancelling()

        states = [
            await self.main_task.get_state(),
            await self.long_running_task.get_state(),
            await self.complete_task.get_state(),
        ]

        if all(s not in PROGRESS_STATES for s in states):
            await self.set_cancelled()
