import uuid
from typing import Type

import pytest
from pydantic import BaseModel

from integration.workflow.mixins import SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin
from tests.integration.conftest import setup
from workflow.core.component import Component, TOutput, TInput
from workflow.core.flow import ComponentFlowHandler
from workflow.core.graph import GraphExecutionNode
from workflow.core.states import ControlStates


class TestComponent(SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin):
    @pytest.mark.asyncio
    async def test_should_successfully_run_component(self, setup, prefork_celery_worker):
        class IO(BaseModel):
            a: int = 10

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[TInput]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return TestComponent

            @property
            def output_parameter_type(self) -> Type[TOutput]:
                return IO

        # arrange
        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        graph = GraphExecutionNode(experiment_id=experiment_id)

        # act
        await graph.schedule(components=[component])
        await graph.execute()
        await self.sync_until_terminal(graph=graph)

        # assert
        assert graph.get_component_node(component_id=component.id).get_state() == ControlStates.SUCCESS
