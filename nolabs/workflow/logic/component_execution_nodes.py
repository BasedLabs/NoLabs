import uuid

from nolabs.infrastructure.log import get_scheduler_logger
from nolabs.infrastructure.red import use_redis_pipe
from nolabs.workflow.logic.graph import CeleryExecutionNode
from nolabs.workflow.logic.states import ControlStates, ERROR_STATES
from nolabs.workflow.logic.graph import ExecutionNode
from nolabs.workflow.logic.tasks import _component_main_task
from nolabs.workflow.logic.job_execution_nodes import JobExecutionNode
from nolabs.workflow.logic.states import TERMINAL_STATES

logger = get_scheduler_logger()


class ComponentMainTaskExecutionNode(CeleryExecutionNode):
    def __init__(self, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{component_id}:main_task")

    async def can_execute(self) -> bool:
        return await self.get_state() == ControlStates.SCHEDULED

    async def execute(self):
        input_data = await self.get_input()
        async with use_redis_pipe() as pipe:
            celery_task_id = await self._prepare_for_execute()
            _component_main_task.apply_async(
                kwargs={"experiment_id": input_data["experiment_id"],
                        "component_id": input_data["component_id"]},
                task_id=celery_task_id,
                retry=False)
            await self.set_state(ControlStates.STARTED)
            await self.set_output(output={})
            await self.set_error(error="")

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
        async with use_redis_pipe():
            input_data = await self.get_input()
            celery_task_id = await self._prepare_for_execute()
            _component_main_task.apply_async(
                kwargs={"experiment_id": input_data["experiment_id"],
                        "component_id": input_data["component_id"],
                        "job_ids": input_data["job_ids"]},
                task_id=celery_task_id,
                retry=False)
            await self.set_state(ControlStates.STARTED)
            await self.set_output(output={})
            await self.set_error(error="")

    async def schedule(self, experiment_id: uuid.UUID, component_id: uuid.UUID):
        async with use_redis_pipe():
            await self.set_input(execution_input={
                'experiment_id': experiment_id,
                'component_id': component_id
            })
            await self.set_state(state=ControlStates.SCHEDULED)


class ComponentExecutionNode(ExecutionNode):
    def __init__(self, component_id: uuid.UUID):
        super().__init__(id=f"execution_node:{component_id}")
        self.component_id = component_id

    async def synchronize(self):
        extra = {
            "component_id": self.component_id
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
        await complete_task.synchronize()

        if state in [ControlStates.SCHEDULED, ControlStates.STARTED]:
            proceed = await self._process_main_task(main_task=main_task)
            if not proceed:
                return
            proceed = await self._process_jobs(main_task=main_task)
            if not proceed:
                return
            await self._process_complete(complete_task=complete_task)

    async def execute(self, **kwargs):
        main_task = ComponentMainTaskExecutionNode(self.component_id)

        if await main_task.can_execute():
            async with use_redis_pipe():
                await self.set_state(state=ControlStates.STARTED)
                await main_task.execute()
                return

        main_task_outputs = main_task.get_output()
        jobs = []
        if main_task_outputs and isinstance(main_task_outputs, list):
            max_batch = 4
            for job_id in main_task_outputs:
                job = JobExecutionNode(job_id=job_id)
                jobs.append((job, await job.get_state()))
            if len([j for j in jobs if j[-1] == ControlStates.STARTED]) >= max_batch:
                return
            executed = False
            async with use_redis_pipe():
                for job, state in jobs:
                    if state == ControlStates.SCHEDULED:
                        await job.execute()
                        executed = True
            if executed:
                return

        if not jobs or all(j for j in jobs if j[-1] in TERMINAL_STATES) and any(
                j for j in jobs if j[-1] == ControlStates.STARTED):
            complete_task = ComponentCompleteTaskExecutionNode(self.component_id)
            await complete_task.execute()

    async def can_execute(self) -> bool:
        pass

    async def _process_main_task(self, main_task: ExecutionNode) -> bool:
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

    async def _process_jobs(self, main_task: ExecutionNode) -> bool:
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

    async def _process_complete(self, complete_task: ExecutionNode) -> bool:
        """
        :returns: move next?
        """
        complete_task_state = await complete_task.get_state()

        if complete_task_state == ControlStates.UNKNOWN:
            inputs = await self.get_input()
            await complete_task.schedule(experiment_id=inputs["experiment_id"], component_id=inputs["component_id"])

        if complete_task_state in TERMINAL_STATES:
            complete_task_output = await complete_task.get_output()
            async with use_redis_pipe():
                await self.set_state(state=complete_task_state)
                await self.set_output(output=complete_task_output)
                if complete_task_state == ControlStates.SUCCESS:
                    return True
                if complete_task_state in ERROR_STATES:
                    error = await complete_task.get_error()
                    await self.set_error(error=error)
        return False