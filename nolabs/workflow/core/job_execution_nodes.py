import uuid
from typing import Optional, Dict, Any

from infrastructure.celery_app_factory import get_celery_app
from infrastructure.log import get_worker_logger
from infrastructure.redis_client_factory import use_redis_pipe
from nolabs.workflow.core.node import CeleryExecutionNode, ExecutionNode
from nolabs.workflow.core.states import ControlStates
from nolabs.workflow.core.celery_tasks import _job_main_task, _complete_job_task
from nolabs.workflow.core.states import ERROR_STATES, TERMINAL_STATES
from workflow.core.socketio_events_emitter import emit_start_job_event, emit_finish_job_event

logger = get_worker_logger()
celery = get_celery_app()


class JobMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self,
                 experiment_id: uuid.UUID,
                 component_id: uuid.UUID,
                 job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}:main_task")
        self.experiment_id = experiment_id
        self.component_id = component_id
        self.job_id = job_id

    async def execute(self):
        async with use_redis_pipe():
            celery_task_id = await self._prepare_for_execute()
            _job_main_task.apply_async(
                kwargs={"experiment_id": self.experiment_id,
                        "component_id": self.component_id,
                        "job_id": self.job_id},
                task_id=celery_task_id,
                retry=False)
            await self.set_state(ControlStates.STARTED)

    async def schedule(self):
        async with use_redis_pipe():
            await self.set_input(execution_input={
                'experiment_id': self.experiment_id,
                'component_id': self.component_id,
                'job_id': self.job_id
            })
            await self.set_state(state=ControlStates.SCHEDULED)
            await self.set_output(output={})
            await self.set_message(message="")


class JobLongRunningTaskExecutionNode(CeleryExecutionNode):
    def __init__(self,
                 experiment_id: uuid.UUID,
                 component_id: uuid.UUID,
                 job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}:long_running_task")

    async def execute(self):
        async with use_redis_pipe():
            input_data = await self.get_input()
            celery_task_name = input_data['celery_task_name']
            kwargs = input_data['kwargs']
            task_id = await self._prepare_for_execute()
            celery.send_task(name=celery_task_name, task_id=task_id, retry=False, kwargs=kwargs)
            await self.set_state(ControlStates.STARTED)
            await self.set_output(output={})
            await self.set_message(message="")

    async def schedule(self, celery_task_name: str, arguments: Optional[Dict[str, Any]] = None):
        data = {
            'celery_task_name': celery_task_name,
            'kwargs': arguments or {}
        }

        async with use_redis_pipe():
            await self.set_input(execution_input=data)
            await self.set_state(state=ControlStates.SCHEDULED)
            await self.set_message(message="")


class JobCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}:complete_task")
        self.job_id = job_id
        self.component_id = component_id
        self.experiment_id = experiment_id

    async def execute(self):
        async with use_redis_pipe():
            input_data = await self.get_input()
            celery_task_id = await self._prepare_for_execute()
            _complete_job_task.apply_async(
                kwargs={"experiment_id": self.experiment_id,
                        "component_id": self.component_id,
                        "job_id": self.job_id,
                        "long_running_output": input_data['long_running_output']
                        },
                task_id=celery_task_id,
                retry=False)
            await self.set_state(ControlStates.STARTED)

    async def schedule(self, long_running_output: Optional[Dict[str, Any]] = None):
        async with use_redis_pipe():
            await self.set_input(execution_input={
                'experiment_id': self.experiment_id,
                'component_id': self.component_id,
                'job_id': self.job_id,
                'long_running_output': long_running_output or {}
            })
            await self.set_state(state=ControlStates.SCHEDULED)
            await self.set_output(output={})
            await self.set_message(message="")


class JobExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}:main_task")
        self.experiment_id = experiment_id
        self.component_id = component_id
        self.job_id = job_id

        self.main_task = JobMainTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id
        )

        self.long_running_task = JobLongRunningTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id
        )

        self.complete_task = JobCompleteTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=self.job_id
        )

    async def sync_started(self):
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
        main_task = JobMainTaskExecutionNode(experiment_id=self.experiment_id,
                                             component_id=self.component_id,
                                             job_id=self.job_id)
        await main_task.sync_started()

        if await main_task.can_schedule():
            await main_task.schedule()
            return

        if await main_task.can_execute():
            await main_task.execute()
            emit_start_job_event(experiment_id=self.experiment_id,
                                 component_id=self.component_id,
                                 job_id=self.job_id)
            return

        main_task_state = await main_task.get_state()
        if main_task_state in TERMINAL_STATES:
            await self.set_state(main_task_state)
            await self.set_message(await main_task.get_message())
            return

    async def sync_long_running_task(self):
        main_task_state = await self.main_task.get_state()

        long_running_task = JobLongRunningTaskExecutionNode(experiment_id=self.experiment_id,
                                                            component_id=self.component_id,
                                                            job_id=self.job_id)
        await long_running_task.sync_started()

        if main_task_state in TERMINAL_STATES and await long_running_task.get_state() == ControlStates.UNKNOWN:
            await self.set_state(main_task_state)
            return

        if await long_running_task.can_execute():
            await long_running_task.execute()
            return

    async def sync_complete_task(self):
        await self.complete_task.sync_started()

        long_running_task = JobLongRunningTaskExecutionNode(experiment_id=self.experiment_id,
                                                            component_id=self.component_id,
                                                            job_id=self.job_id)
        long_running_task_state = await long_running_task.get_state()

        if long_running_task_state == ControlStates.SUCCESS and await self.complete_task.can_schedule():
            await self.complete_task.schedule(long_running_output=await long_running_task.get_output())
            return

        if await self.complete_task.can_execute():
            await self.complete_task.execute()
            return

    async def set_final_state(self):
        success = True

        for task in [self.main_task, self.long_running_task, self.complete_task]:
            state = await task.get_state()
            if state in ERROR_STATES:
                success = False
                await self.set_state(state)
                await self.set_message(await task.get_message())
                emit_finish_job_event(experiment_id=self.experiment_id,
                                      component_id=self.component_id,
                                      job_id=self.job_id)
                return

        if success:
            await self.set_state(ControlStates.SUCCESS)
            await self.set_message("")
            emit_finish_job_event(experiment_id=self.experiment_id,
                                  component_id=self.component_id,
                                  job_id=self.job_id)

    async def execute(self):
        await self.set_state(state=ControlStates.STARTED)

    async def schedule(self,
                       experiment_id: uuid.UUID,
                       component_id: uuid.UUID,
                       job_id: uuid.UUID):
        await self.set_input(execution_input={
            'experiment_id': experiment_id,
            'component_id': component_id,
            'job_id': job_id
        })
        await self.set_state(state=ControlStates.SCHEDULED)
