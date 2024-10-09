import uuid
from typing import Optional, Dict, Any

from domain.exceptions import NoLabsException, ErrorCodes
from infrastructure.cel import cel
from infrastructure.log import get_scheduler_logger
from infrastructure.red import use_redis_pipe
from workflow.logic.control import _job_main_task, _complete_job_task
from workflow.logic.graph import CeleryExecutionNode, ExecutionNode
from workflow.logic.states import ControlStates, TERMINAL_STATES, ERROR_STATES

logger = get_scheduler_logger()


class JobMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{job_id}:main_task")

    async def can_execute(self) -> bool:
        return await self.get_state() == ControlStates.SCHEDULED

    async def execute(self):
        input_data = await self.get_input()
        celery_task_id = await self._prepare_for_execute()
        _job_main_task.apply_async(
            kwargs={"experiment_id": input_data["experiment_id"],
                    "component_id": input_data["component_id"],
                    "job_id": input_data["job_id"]},
            task_id=celery_task_id,
            retry=False)

    async def schedule(self, experiment_id: uuid.UUID, component_id: uuid.UUID, job_id: uuid.UUID):
        async with use_redis_pipe():
            await self.set_input(execution_input={
                'experiment_id': experiment_id,
                'component_id': component_id,
                'job_id': job_id
            })
            await self.set_state(state=ControlStates.SCHEDULED)


class JobLongRunningTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{job_id}:long_running_task")

    async def can_execute(self) -> bool:
        return await self.get_state() == ControlStates.SCHEDULED

    async def execute(self):
        input_data = await self.get_input()
        celery_task_name = input_data['celery_task_name']
        kwargs = input_data['kwargs']
        task_id = await self._prepare_for_execute()
        cel.app.send_task(name=celery_task_name, task_id=task_id, retry=False, kwargs=kwargs)

    async def schedule(self, celery_task_name: str, arguments: Optional[Dict[str, Any]] = None):
        data = {
            'celery_task_name': celery_task_name,
            'kwargs': arguments or {}
        }

        async with use_redis_pipe():
            await self.set_input(execution_input=data)
            await self.set_state(state=ControlStates.SCHEDULED)


class JobCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{job_id}:complete_task")

    async def can_execute(self) -> bool:
        return await self.get_state() == ControlStates.SCHEDULED

    async def execute(self):
        input_data = await self.get_input()
        celery_task_id = await self._prepare_for_execute()
        _complete_job_task.apply_async(
            kwargs={"experiment_id": input_data["experiment_id"],
                    "component_id": input_data["component_id"],
                    "job_id": input_data["job_id"],
                    "long_running_output": input_data["long_running_output"]
                    },
            task_id=celery_task_id,
            retry=False)

    async def schedule(self,
                       experiment_id: uuid.UUID,
                       component_id: uuid.UUID,
                       job_id: uuid.UUID,
                       long_running_output: Optional[Dict[str, Any]] = None):
        if not experiment_id or not component_id or not job_id:
            raise NoLabsException(ErrorCodes.cannot_schedule_execution_node,
                                  message="kwargs fields: experiment_id, component_id, job_id, long_running_output")

        async with use_redis_pipe():
            await self.set_input(execution_input={
                'experiment_id': experiment_id,
                'component_id': component_id,
                'job_id': job_id,
                'long_running_output': long_running_output or {}
            })
            await self.set_state(state=ControlStates.SCHEDULED)


class JobExecutionNode(ExecutionNode):
    def __init__(self, job_id: uuid.UUID):
        super().__init__(id=f"execution_node:{job_id}")
        self.job_id = job_id

    async def synchronize(self):
        async def process_main_task() -> bool:
            """
            :returns: move next?
            """
            main_task_state = await main_task.get_state()

            if main_task_state == ControlStates.UNKNOWN:
                inputs = await self.get_input()
                await main_task.schedule(experiment_id=inputs["experiment_id"],
                                         component_id=inputs["component_id"],
                                         job_id=inputs["job_id"])
                return False

            if main_task_state == ControlStates.SUCCESS:
                long_running_task_state = await long_running_task.get_state()
                if long_running_task_state == ControlStates.UNKNOWN:
                    # We did not schedule long-running task
                    await self.set_state(state=ControlStates.SUCCESS)
                    return False
                return True
            if main_task_state in ERROR_STATES:
                async with use_redis_pipe():
                    error = await main_task.get_error()
                    await self.set_state(state=main_task_state)
                    await self.set_error(error=error)
                    await self.set_output(output={})
            return False

        async def process_long_running() -> bool:
            """
            This function is intended to start after process_main_task.
            Assuming that main_task_state is OK
            :returns: move next?
            """
            long_running_task_state = await long_running_task.get_state()

            if long_running_task_state == ControlStates.SUCCESS:
                return True

            if long_running_task_state in ERROR_STATES:
                async with use_redis_pipe():
                    error = await long_running_task.get_error()
                    await self.set_state(state=long_running_task_state)
                    await self.set_error(error=error)
                    await self.set_output(output={})
            return False

        async def process_complete() -> bool:
            """
            :returns: move next?
            """
            complete_task_state = await complete_task.get_state()

            if complete_task_state in TERMINAL_STATES:
                complete_task_output = await complete_task.get_output()
                async with use_redis_pipe():
                    await self.set_state(state=complete_task_state)
                    await self.set_output(output=complete_task_output)
                    if complete_task_state == ControlStates.SUCCESS:
                        return True
                    if complete_task_state in ERROR_STATES:
                        error = await long_running_task.get_error()
                        await self.set_error(error=error)
            return False

        extra = {
            "job_id": self.job_id
        }
        state = await self.get_state()

        if state in (TERMINAL_STATES + [ControlStates.UNKNOWN]):
            return

        main_task = JobMainTaskExecutionNode(self.job_id)
        long_running_task = JobLongRunningTaskExecutionNode(self.job_id)
        complete_task = JobCompleteTaskExecutionNode(self.job_id)

        if state == ControlStates.SCHEDULED:
            inputs = await self.get_input()

            if await main_task.get_state() == ControlStates.UNKNOWN:
                await main_task.schedule(experiment_id=inputs["experiment_id"],
                                         component_id=inputs["component_id"],
                                         job_id=inputs["job_id"])
                return

            logger.critical("Job execution node is scheduled, but main task is not in unknown state", extra=extra)
            return

        await main_task.synchronize()
        await long_running_task.synchronize()
        await complete_task.synchronize()

        if state in [ControlStates.SCHEDULED, ControlStates.STARTED]:
            proceed = await process_main_task()
            if not proceed:
                return
            proceed = await process_long_running()
            if not proceed:
                return
            await process_complete()

    async def can_execute(self) -> bool:
        t1 = JobMainTaskExecutionNode(self.job_id)
        t2 = JobLongRunningTaskExecutionNode(self.job_id)
        t3 = JobCompleteTaskExecutionNode(self.job_id)

        return await t1.can_execute() or await t2.can_execute() or t3.can_execute()

    async def execute(self, **kwargs):
        tasks = [
            JobMainTaskExecutionNode(self.job_id),
            JobLongRunningTaskExecutionNode(self.job_id),
            JobCompleteTaskExecutionNode(self.job_id)
        ]

        for task in tasks:
            if await task.get_state() == ControlStates.SCHEDULED:
                await task.execute()
                return

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
