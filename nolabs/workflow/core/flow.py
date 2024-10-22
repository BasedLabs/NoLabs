import asyncio
import threading
import uuid
from typing import Generic, List, Optional, Dict, Any, Union

from pydantic import BaseModel

from nolabs.workflow.core.component import TInput, TOutput
from nolabs.workflow.core.job_execution_nodes import JobLongRunningTaskExecutionNode


class ComponentFlowHandler(Generic[TInput, TOutput]):
    """
    Represents component flow and available client handlers
    """

    def __init__(self,
                 component_id: uuid.UUID,
                 experiment_id: uuid.UUID):
        self.component_id = component_id
        self.experiment_id = experiment_id

    async def on_component_task(self, inp: TInput) -> List[uuid.UUID]:
        return []

    async def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
        return None

    async def on_job_task(self, job_id: uuid.UUID):
        pass

    async def on_job_completion(self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]):
        pass

    async def schedule_long_running(self, job_id: uuid.UUID, celery_task_name: Any, input: Optional[Union[BaseModel, Dict[str, Any]]] = None):
        node = JobLongRunningTaskExecutionNode(experiment_id=self.experiment_id, component_id=self.component_id, job_id=job_id)
        can_schedule = await node.can_schedule()
        if can_schedule:
            arguments = {}

            if isinstance(input, dict):
                arguments = input
            elif input:
                arguments = input.model_dump()

            await node.schedule(celery_task_name=celery_task_name, arguments=arguments)
