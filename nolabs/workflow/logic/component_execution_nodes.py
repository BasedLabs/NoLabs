import uuid

from nolabs.infrastructure.log import get_scheduler_logger
from nolabs.infrastructure.red import use_redis_pipe
from nolabs.workflow.logic.graph import CeleryExecutionNode
from nolabs.workflow.logic.states import ControlStates, ERROR_STATES
from nolabs.workflow.logic.graph import _component_main_task
from workflow.logic.graph import ExecutionNode
from workflow.logic.job_execution_nodes import JobExecutionNode
from workflow.logic.states import TERMINAL_STATES

logger = get_scheduler_logger()


class ComponentMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{component_id}:main_task")

    async def can_execute(self) -> bool:
        return await self.get_state() == ControlStates.SCHEDULED

    async def execute(self):
        input_data = await self.get_input()
        celery_task_id = await self._prepare_for_execute()
        _component_main_task.apply_async(
            kwargs={"experiment_id": input_data["experiment_id"],
                    "component_id": input_data["component_id"]},
            task_id=celery_task_id,
            retry=False)

    async def schedule(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        async with use_redis_pipe():
            await self.set_input(execution_input={
                'experiment_id': experiment_id,
                'component_id': component_id
            })
            await self.set_state(state=ControlStates.SCHEDULED)


class ComponentCompleteTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{component_id}:complete_task")

    async def can_execute(self) -> bool:
        return await self.get_state() == ControlStates.SCHEDULED

    async def execute(self):
        input_data = await self.get_input()
        celery_task_id = await self._prepare_for_execute()
        _component_main_task.apply_async(
            kwargs={"experiment_id": input_data["experiment_id"],
                    "component_id": input_data["component_id"],
                    "job_ids": input_data["job_ids"]},
            task_id=celery_task_id,
            retry=False)

    async def schedule(self, **kwargs):
        experiment_id = kwargs.get("experiment_id")
        component_id = kwargs.get("component_id")

        async with use_redis_pipe():
            await self.set_input(execution_input={
                'experiment_id': experiment_id,
                'component_id': component_id
            })
            await self.set_state(state=ControlStates.SCHEDULED)


class ComponentTaskExecutionNode(ExecutionNode):
    def __init__(self, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{component_id}")
        self.component_id = component_id

    async def synchronize(self):
        async def process_main_task() -> bool:
            """
            :returns: propagate
            """
            main_task_state = await main_task.get_state()

            if main_task_state == ControlStates.UNKNOWN:
                inputs = await self.get_input()
                await main_task.schedule(experiment_id=inputs["experiment_id"],
                                         component_id=inputs["component_id"])
                return False

            if main_task_state == ControlStates.SUCCESS:
                return True

            if main_task_state in ERROR_STATES:
                async with use_redis_pipe():
                    error = await main_task.get_error()
                    await self.set_state(state=main_task_state)
                    await self.set_error(error=error)
                    await self.set_output(output={})
            return False

        async def process_jobs() -> bool:
            """
            This function is intended to start after process_main_task.
            Assuming that main_task_state is OK
            :returns: propagate
            """
            main_task_outputs = main_task.get_output()
            if main_task_outputs and isinstance(main_task_outputs, list):
                component_inputs = await self.get_input()
                job_ids = main_task_outputs
                states = []
                for job_id in job_ids:
                    job = JobExecutionNode(job_id=job_id)
                    await job.synchronize()
                    job_state = await job.get_state()
                    states.append(job_state)
                    if job_state == ControlStates.UNKNOWN:
                        await job.schedule(experiment_id=component_inputs["experiment_id"],
                                           component_id=component_inputs["component_id"],
                                           job_id=job_id
                                           )

                if any(s in [ControlStates.SCHEDULED, ControlStates.STARTED] for s in states):
                    return False

                if any(s == ControlStates.SUCCESS for s in states):
                    return True

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

        main_task = ComponentMainTaskExecutionNode(self.component_id)
        complete_task = ComponentCompleteTaskExecutionNode(self.component_id)

        if state == ControlStates.SCHEDULED:
            inputs = await self.get_input()

            if await main_task.get_state() == ControlStates.UNKNOWN:
                await main_task.schedule(experiment_id=inputs["experiment_id"],
                                         component_id=inputs["component_id"])
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

    async def execute(self, **kwargs):
        pass

    async def can_execute(self) -> bool:
        pass
