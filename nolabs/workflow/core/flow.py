import uuid
from typing import Any, Dict, Generic, List, Optional, Union

from pydantic import BaseModel

from nolabs.workflow.core.component import TInput, TOutput
from nolabs.workflow.core.job_execution_nodes import JobLongRunningTaskExecutionNode, JobMainTaskExecutionNode, \
    JobExecutionNode


class ComponentFlowHandler(Generic[TInput, TOutput]):
    """
    Represents component flow and available client handlers
    """

    def __init__(self, component_id: uuid.UUID, experiment_id: uuid.UUID):
        self.component_id = component_id
        self.experiment_id = experiment_id

    async def on_start(self, inp: TInput) -> List[uuid.UUID]:
        return []

    async def on_finish(
            self, inp: TInput, job_ids: List[uuid.UUID]
    ) -> Optional[TOutput]:
        return None

    async def on_job_start(self, job_id: uuid.UUID):
        pass

    async def on_job_finish(
            self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
    ):
        pass

    async def cancel_job(self, job_id: uuid.UUID, reason: Optional[str] = None):
        node = JobExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=job_id,
        )
        if await node.can_cancel():
            await node.cancel(message=reason)

    async def schedule(
            self,
            job_id: uuid.UUID,
            celery_task_name: str,
            celery_queue: str,
            input: Optional[Union[BaseModel, Dict[str, Any]]] = None
    ):
        node = JobLongRunningTaskExecutionNode(
            experiment_id=self.experiment_id,
            component_id=self.component_id,
            job_id=job_id,
        )
        can_schedule = await node.can_schedule()
        if can_schedule:
            arguments = {}

            if isinstance(input, dict):
                arguments = input
            elif input:
                arguments = input.model_dump()

            await node.schedule(
                celery_task_name=celery_task_name,
                arguments=arguments,
                celery_queue=celery_queue,
            )
