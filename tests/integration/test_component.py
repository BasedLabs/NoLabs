import uuid
from typing import Type, List, Optional

from pydantic import BaseModel

from integration.mixins import SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin
from integration.setup import GlobalSetup
from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler
from nolabs.workflow.core.graph import GraphExecutionNode
from nolabs.workflow.core.states import ControlStates


class TestComponent(GlobalSetup,
               SeedComponentsMixin,
               SeedExperimentMixin,
               GraphTestMixin):
    async def test_should_successfully_run_component(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            ...

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[TInput]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return FlowHandler

            @property
            def output_parameter_type(self) -> Type[TOutput]:
                return IO

        # arrange
        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        graph = GraphExecutionNode(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.schedule(components=[component])
        await graph.start()
        await self.sync_until_terminal(graph=graph)

        # assert
        self.assertEqual(await graph.get_component_node(component_id=component.id).get_state(), ControlStates.SUCCESS)


    async def test_should_fail_component_on_main_task_failure(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            def on_component_task(self, inp: TInput) -> List[uuid.UUID]:
                raise ValueError("Hello")

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[TInput]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return FlowHandler

            @property
            def output_parameter_type(self) -> Type[TOutput]:
                return IO

        # arrange
        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        graph = GraphExecutionNode(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.schedule(components=[component])
        await graph.start()
        await self.sync_until_terminal(graph=graph)

        # assert
        self.assertEqual(await graph.get_component_node(component_id=component.id).get_state(), ControlStates.FAILURE)
        self.assertEqual(await graph.get_component_node(component_id=component.id).get_message(), "Hello")

    async def test_should_fail_component_on_complete_task_failure(self):
        # arrange

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
                raise ValueError("Hello")

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[TInput]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return FlowHandler

            @property
            def output_parameter_type(self) -> Type[TOutput]:
                return IO

        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        graph = GraphExecutionNode(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.schedule(components=[component])
        await graph.start()
        await self.sync_until_terminal(graph=graph)

        # assert
        self.assertEqual(await graph.get_component_node(component_id=component.id).get_state(), ControlStates.FAILURE)
        self.assertEqual(await graph.get_component_node(component_id=component.id).get_message(), "Hello")


