import uuid
from typing import List, Optional, Dict

from pydantic import BaseModel

from nolabs.infrastructure.redis_client_factory import get_redis_pipe
from nolabs.workflow.core.component import Component
from nolabs.workflow.core.component_execution_nodes import ComponentExecutionNode
from nolabs.workflow.core.job_execution_nodes import JobExecutionNode
from nolabs.workflow.core.node import ExecutionNode
from nolabs.workflow.core.states import ControlStates
from nolabs.workflow.core.states import TERMINAL_STATES


class GraphData(BaseModel):
    graph: Dict[uuid.UUID, List[uuid.UUID]]
    schedule: List[uuid.UUID]



class GraphExecutionNode(ExecutionNode):
    def __init__(self, experiment_id: uuid.UUID):
        super().__init__(f"execution_node:{experiment_id}")
        self.experiment_id = experiment_id

    async def sync_started(self):
        model = GraphData(**await self.get_input())
        graph = model.graph
        schedule = model.schedule

        all_completed = True

        for component_id in graph.keys():
            node = ComponentExecutionNode(experiment_id=self.experiment_id, component_id=component_id)

            previous_component_ids = graph[component_id]

            if component_id in schedule:
                if await node.can_schedule(previous_component_ids=previous_component_ids):
                    await node.schedule()

            if await node.can_start():
                await node.start()

            await node.sync_started()
            node_state = await node.get_state()
            # By checking for Unknown state we skipped not scheduled components
            if node_state not in TERMINAL_STATES and node_state != ControlStates.UNKNOWN:
                all_completed = False

        if all_completed:
            await self.set_state(ControlStates.SUCCESS)

    async def start(self, **kwargs):
        await self.set_state(ControlStates.STARTED)

    async def schedule(self, components: List[Component], schedule: Optional[List[Component]] = None):
        graph = {}

        if not schedule:
            schedule = components

        for component in components:
            graph[component.id] = []
            for component_2 in components:
                if component_2.id == component.id:
                    continue
                if component_2.id in component.previous_component_ids:
                    graph[component.id].append(component_2.id)

        pipe = get_redis_pipe()
        model = GraphData(
            graph=graph,
            schedule=[c.id for c in schedule]
        )
        await self.set_input(execution_input=model, pipe=pipe)
        await self.set_state(state=ControlStates.SCHEDULED, pipe=pipe)
        await pipe.execute()

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
