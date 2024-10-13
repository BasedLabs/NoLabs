import json
import uuid
from typing import List, Set, Dict, Optional

from infrastructure.log import get_worker_logger
from nolabs.workflow.core.component import Component
from nolabs.workflow.core.states import ControlStates
from nolabs.workflow.core.component_execution_nodes import ComponentExecutionNode
from nolabs.workflow.core.node import ExecutionNode
from infrastructure.redis_client_factory import use_redis_pipe, get_redis_client
from workflow.core.job_execution_nodes import JobExecutionNode

logger = get_worker_logger()
redis_client = get_redis_client()

class GraphExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID):
        super().__init__(f"execution_node:{experiment_id}")
        self.experiment_id = experiment_id

    async def sync_started(self):
        graph = await self.get_input()
        for component_id in graph.keys():
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=uuid.UUID(component_id))
            await node.sync_started()

    async def execute(self, **kwargs):
        graph = await self.get_input()
        for component_id in graph.keys():
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=uuid.UUID(component_id))
            await node.execute()

    async def schedule(self, components: List[Component]):
        graph = {}

        for component in components:
            graph[str(component.id)] = set()
            for component_2 in components:
                if component_2.id == component.id:
                    continue
                if component_2.id in component.previous_component_ids:
                    graph[str(component.id)].add(component_2.id)

        async with use_redis_pipe():
            await self.set_input(execution_input=graph)
            await self.set_state(state=ControlStates.SCHEDULED)

    def get_component_node(self, component_id: uuid.UUID) -> ComponentExecutionNode:
        return ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component_id)

    def get_job_node(self, component_id: uuid.UUID, job_id: uuid.UUID) -> JobExecutionNode:
        return JobExecutionNode(experiment_id=self.experiment_id, component_id=component_id, job_id=job_id)

    async def can_schedule_job(self, component_id: uuid.UUID, job_id: uuid.UUID) -> bool:
        node = JobExecutionNode(experiment_id=self.experiment_id, component_id=component_id, job_id=job_id)
        return await node.can_schedule()

    async def schedule_job(self, component_id: uuid.UUID, job_id: uuid.UUID):
        node = JobExecutionNode(experiment_id=self.experiment_id, component_id=component_id, job_id=job_id)
        await node.schedule(experiment_id=self.experiment_id, component_id=component_id, job_id=job_id)
