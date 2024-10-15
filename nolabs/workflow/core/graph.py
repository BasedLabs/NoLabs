import uuid
from typing import List

from infrastructure.redis_client_factory import use_redis_pipe
from nolabs.workflow.core.component import Component
from nolabs.workflow.core.component_execution_nodes import ComponentExecutionNode
from nolabs.workflow.core.node import ExecutionNode
from nolabs.workflow.core.states import ControlStates
from nolabs.workflow.core.job_execution_nodes import JobExecutionNode
from nolabs.workflow.core.states import TERMINAL_STATES


class GraphExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID):
        super().__init__(f"execution_node:{experiment_id}")
        self.experiment_id = experiment_id

    async def sync_started(self):
        graph = await self.get_input()

        all_completed = True

        for component_id in graph.keys():
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=uuid.UUID(component_id))
            await node.sync_started()
            if await node.get_state() not in TERMINAL_STATES:
                all_completed = False

        if all_completed:
            await self.set_state(ControlStates.SUCCESS)

    async def start(self, **kwargs):
        await self.set_state(ControlStates.STARTED)
        graph = await self.get_input()
        for component_id in graph.keys():
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=uuid.UUID(component_id))
            await node.start()

    async def schedule(self, components: List[Component]):
        graph = {}

        for component in components:
            graph[str(component.id)] = []
            for component_2 in components:
                if component_2.id == component.id:
                    continue
                if component_2.id in component.previous_component_ids:
                    graph[str(component.id)].append(component_2.id)

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
