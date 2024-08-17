from typing import List
from uuid import UUID

from nolabs.application.workflow.api import Component
from nolabs.application.workflow.dag_generator import generate_workflow_dag
from nolabs.infrastructure.environment import Environment
from nolabs.infrastructure.settings import Settings


class DagExecutor:
    _settings: Settings

    def __init__(self):
        self._settings = Settings.load()

    async def execute(self, workflow_id: UUID, experiment_id: UUID, components: List[Component]):
        dag = generate_workflow_dag(workflow_id=workflow_id, experiment_id=experiment_id, components=components)

        if self._settings.get_environment() == Environment.LOCAL:
            dag.test()

