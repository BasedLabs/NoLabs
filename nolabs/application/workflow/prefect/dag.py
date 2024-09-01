import uuid
from typing import List, Optional, Dict, Any
from uuid import UUID

from prefect import flow

from nolabs.application.workflow.workflow import Component
from nolabs.infrastructure.environment import Environment
from nolabs.infrastructure.settings import Settings


class PrefectDagExecutor:
    _settings: Settings

    def __init__(self):
        self._settings = Settings.load()

    async def execute(self, workflow_id: UUID, components: List[Component], extra: Optional[Dict[str, Any]] = None):
        dag = generate_workflow_dag(workflow_id=workflow_id, components=components, extra=extra)

        if self._settings.get_environment() == Environment.LOCAL:
            await dag()


def generate_workflow_dag(workflow_id: uuid.UUID, components: List[Component], extra: Optional[Dict[str, Any]] = None):
    def component_flow_factory(component: Component):
        @flow(name=str(component.name), flow_run_name=f'{component.name}-{str(component.id)}')
        async def component_flow():
            child_flows = []

            for child_component_id in component.previous_component_ids:
                child_component = [c for c in components if c.id == child_component_id][0]
                child_flow = component_flow_factory(component=child_component)
                child_flow_result = await child_flow()
                child_flows.append(child_flow_result)

            component_task_type = component.component_task_type(
                component_id=component.id,
                component_name=component.name,
                workflow_id=workflow_id,
                extra=extra
            )

            await component_task_type.execute_task()

        return component_flow

    @flow(name=str(workflow_id), flow_run_name=f'workflow-{str(workflow_id)}')
    async def workflow():
        for component in components:
            await component_flow_factory(component=component)()

    return workflow
