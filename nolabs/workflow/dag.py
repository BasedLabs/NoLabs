import asyncio
import uuid
from collections import defaultdict, deque
from typing import Dict, List, Set
from uuid import UUID

from infrastructure.settings import settings
from prefect import State, flow
from prefect.client.schemas import FlowRun
from workflow.component import Component
from workflow.data import WorkflowData

from nolabs.infrastructure.environment import Environment


class PrefectDagExecutor:
    async def execute(
        self, workflow_id: UUID, experiment_id: uuid.UUID, components: List[Component]
    ):
        dag = generate_workflow_dag(
            workflow_id=workflow_id, components=components, experiment_id=experiment_id
        )

        if settings.get_environment() == Environment.LOCAL:
            await dag(return_state=True)


def generate_workflow_dag(
    workflow_id: uuid.UUID, experiment_id: uuid.UUID, components: List[Component]
):
    def on_running(_, flow_run: FlowRun, state: State):
        data: WorkflowData = WorkflowData.objects.with_id(workflow_id)
        data.flow_run_id = flow_run.id
        data.save()

    @flow(
        name=str(workflow_id),
        flow_run_name=f"Workflow,{str(workflow_id)}",
        on_running=[on_running],
    )
    async def workflow():
        if not components:
            return

        if len(components) == 1:
            component = components[0]
            component_flow = component.component_flow_type(
                component=component, experiment_id=experiment_id
            )
            await component_flow.execute(return_state=True)
            return

        component_map: Dict[uuid.UUID, Component] = {c.id: c for c in components}

        in_degree: Dict[uuid.UUID, int] = defaultdict(int)
        adj_list: Dict[uuid.UUID, Set[uuid.UUID]] = defaultdict(set)

        for component in components:
            for prev_id in component.previous_component_ids:
                adj_list[prev_id].add(component.id)
                in_degree[component.id] += 1

        queue = deque([c.id for c in components if in_degree[c.id] == 0])

        while queue:
            parallel_group = []

            for _ in range(len(queue)):
                current_id = queue.popleft()
                current_component = component_map[current_id]
                parallel_group.append(current_component)

            group = [
                c.component_flow_type(component=c, experiment_id=experiment_id).execute(
                    return_state=True
                )
                for c in parallel_group
            ]

            await asyncio.gather(*group)

            for component in parallel_group:
                for dependent_id in adj_list[component.id]:
                    in_degree[dependent_id] -= 1

                    if in_degree[dependent_id] == 0:
                        queue.append(dependent_id)

    return workflow
