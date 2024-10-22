import uuid
from typing import Optional, Dict, Any

from nolabs.infrastructure.redis_client_factory import get_redis_pipe
from nolabs.workflow.core.celery_tasks import Tasks
from nolabs.workflow.core.node import CeleryExecutionNode, ExecutionNode
from nolabs.workflow.core.socketio_events_emitter import emit_finish_job_event
from nolabs.workflow.core.states import ControlStates
from nolabs.workflow.core.states import ERROR_STATES, TERMINAL_STATES


class JobMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self,
                 experiment_id: uuid.UUID,
                 component_id: uuid.UUID,
                 job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}:main_task")
        self.experiment_id = experiment_id
        self.component_id = component_id
        self.job_id = job_id

    async def start(self):
        pipe = get_redis_pipe()
        celery_task_id = await self._prepare_for_start(pipe=pipe)
        self.celery.send_task(name=Tasks.job_main_task,
                              kwargs={"experiment_id": self.experiment_id,
                                      "component_id": self.component_id,
                                      "job_id": self.job_id},
                              task_id=celery_task_id,
                              retry=False
                              )
        await self.set_state(ControlStates.STARTED,pipe=pipe)
        await pipe.execute()

    async def schedule(self):
        pipe = get_redis_pipe()
        await self.set_input(execution_input={
            'experiment_id': self.experiment_id,
            'component_id': self.component_id,
            'job_id': self.job_id
        }, pipe=pipe)
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        await self.set_output(output={}, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        await pipe.execute()


class JobLongRunningTaskExecutionNode(CeleryExecutionNode):
    def __init__(self,
                 experiment_id: uuid.UUID,
                 component_id: uuid.UUID,
                 job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}:long_running_task")

    async def start(self):
        input_data = await self.get_input()
        pipe = get_redis_pipe()
        celery_task_name = input_data['celery_task_name']
        kwargs = input_data['kwargs']
        task_id = await self._prepare_for_start(pipe=pipe)
        self.celery.send_task(name=celery_task_name, task_id=task_id, retry=False, kwargs=kwargs)
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        await self.set_output(output={}, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        await pipe.execute()

    async def schedule(self, celery_task_name: str, arguments: Optional[Dict[str, Any]] = None):
        data = {
            'celery_task_name': celery_task_name,
            'kwargs': arguments or {}
        }

        pipe = get_redis_pipe()
        await self.set_input(execution_input=data, pipe=pipe)
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        await pipe.execute()


class JobCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}:complete_task")
        self.job_id = job_id
        self.component_id = component_id
        self.experiment_id = experiment_id

    async def start(self):
        input_data = await self.get_input()
        pipe = get_redis_pipe()
        celery_task_id = await self._prepare_for_start(pipe=pipe)
        self.celery.send_task(name=Tasks.complete_job_task,
                              kwargs={"experiment_id": self.experiment_id,
                                      "component_id": self.component_id,
                                      "job_id": self.job_id,
                                      "long_running_output": input_data['long_running_output']
                                      },
                              task_id=celery_task_id,
                              retry=False
                              )
        await self.set_state(ControlStates.STARTED, pipe=pipe)
        await pipe.execute()

    async def schedule(self, long_running_output: Optional[Dict[str, Any]] = None):
        pipe = get_redis_pipe()
        await self.set_input(execution_input={
            'experiment_id': self.experiment_id,
            'component_id': self.component_id,
            'job_id': self.job_id,
            'long_running_output': long_running_output or {}
        }, pipe=pipe)
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        await self.set_output(output={}, pipe=pipe)
        await self.set_message(message="", pipe=pipe)
        await pipe.execute()


class JobExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{experiment_id}:{component_id}:{job_id}")
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
        if await self.main_task.can_schedule():
            await self.main_task.schedule()

        if await self.main_task.can_start():
            await self.main_task.start()

        await self.main_task.sync_started()

    async def sync_long_running_task(self):
        main_task_state = await self.main_task.get_state()

        # will not run
        if main_task_state in TERMINAL_STATES and await self.long_running_task.get_state() == ControlStates.UNKNOWN:
            return

        if await self.long_running_task.can_start():
            await self.long_running_task.start()
            return

        await self.long_running_task.sync_started()

    async def sync_complete_task(self):
        long_running_task_state = await self.long_running_task.get_state()
        main_task_state = await self.main_task.get_state()

        if main_task_state == ControlStates.SUCCESS and long_running_task_state in [ControlStates.SUCCESS,
                                                                                    ControlStates.UNKNOWN]:
            if await self.complete_task.can_schedule():
                await self.complete_task.schedule(long_running_output=await self.long_running_task.get_output())
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
                emit_finish_job_event(experiment_id=self.experiment_id,
                                      component_id=self.component_id,
                                      job_id=self.job_id)
                break

        if success:
            await self.set_state(ControlStates.SUCCESS)
            await self.set_message("")

        emit_finish_job_event(experiment_id=self.experiment_id,
                              component_id=self.component_id,
                              job_id=self.job_id)

    async def start(self):
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
