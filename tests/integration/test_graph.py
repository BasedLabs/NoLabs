import uuid
from typing import Type, List, Optional

from pydantic import BaseModel

from integration.mixins import SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin
from integration.setup import GlobalSetup
from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler
from nolabs.workflow.core.graph import GraphExecutionNode
from nolabs.workflow.core.states import ControlStates
from nolabs.workflow.core.syncer import Syncer


class TestGraph(GlobalSetup,
               SeedComponentsMixin,
               SeedExperimentMixin,
               GraphTestMixin):
    async def test_should_complete_graph_if_component_failed(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            def on_component_task(self, inp: TInput) -> List[uuid.UUID]:
                raise ValueError('Exception')

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
        syncer = Syncer()
        await syncer.sync_graph(experiment_id=experiment_id, wait=True)

        # assert
        self.assertEqual(await graph.get_state(), ControlStates.SUCCESS)

    async def test_should_run_linear_workflow(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
                return IO(a=15)

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
        component1 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        component2 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        self.seed_mappings(component2, previous_components=[(component1, ['a'], ['a'])])
        graph = GraphExecutionNode(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.schedule(components=[component1, component2])
        syncer = Syncer()
        await syncer.sync_graph(experiment_id=experiment_id, wait=True)

        # assert
        self.assertEqual(await graph.get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component1.id).get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component2.id).get_state(), ControlStates.SUCCESS)

    async def test_should_run_diamond_workflow(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
                return IO(a=15)

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
        component1 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        component2 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        component3 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        component4 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        self.seed_mappings(component2, previous_components=[(component1, ['a'], ['a'])])
        self.seed_mappings(component3, previous_components=[(component1, ['a'], ['a'])])
        self.seed_mappings(component4, previous_components=[(component2, ['a'], ['a']), (component3, ['a'], ['a'])])
        graph = GraphExecutionNode(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.schedule(components=[component1, component2, component3, component4])
        syncer = Syncer()
        await syncer.sync_graph(experiment_id=experiment_id, wait=True)

        # assert
        self.assertEqual(await graph.get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component1.id).get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component2.id).get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component3.id).get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component4.id).get_state(), ControlStates.SUCCESS)

    async def test_should_not_propagate_execution_if_component_failed(self):
        component2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
                if self.component_id == component2_id:
                    raise ValueError("Component failed to propagate")
                return IO(a=15)

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
        component1 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        component2 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent, component_id=component2_id)
        component3 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        component4 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        self.seed_mappings(component2, previous_components=[(component1, ['a'], ['a'])])
        self.seed_mappings(component3, previous_components=[(component1, ['a'], ['a'])])
        self.seed_mappings(component4, previous_components=[(component2, ['a'], ['a']), (component3, ['a'], ['a'])])
        graph = GraphExecutionNode(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.schedule(components=[component1, component2, component3, component4])
        syncer = Syncer()
        await syncer.sync_graph(experiment_id=experiment_id, wait=True)

        # assert
        self.assertEqual(await graph.get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component1.id).get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component2.id).get_state(), ControlStates.FAILURE)
        self.assertEqual(await graph.get_component_node(component_id=component3.id).get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component4.id).get_state(), ControlStates.UNKNOWN)

    async def test_should_not_propagate_execution_if_component_failed(self):
        component2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_completion(self, inp: TInput, job_ids: List[uuid.UUID]) -> Optional[TOutput]:
                if self.component_id == component2_id:
                    raise ValueError("Component failed to propagate")
                return IO(a=15)

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
        component1 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        component2 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent, component_id=component2_id)
        component3 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        component4 = self.seed_component(experiment_id=experiment_id, component_type=MockComponent)
        self.seed_mappings(component2, previous_components=[(component1, ['a'], ['a'])])
        self.seed_mappings(component3, previous_components=[(component1, ['a'], ['a'])])
        self.seed_mappings(component4, previous_components=[(component2, ['a'], ['a']), (component3, ['a'], ['a'])])
        graph = GraphExecutionNode(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.schedule(components=[component1, component2, component3, component4])
        syncer = Syncer()
        await syncer.sync_graph(experiment_id=experiment_id, wait=True)

        # assert
        self.assertEqual(await graph.get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component1.id).get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component2.id).get_state(), ControlStates.FAILURE)
        self.assertEqual(await graph.get_component_node(component_id=component3.id).get_state(), ControlStates.SUCCESS)
        self.assertEqual(await graph.get_component_node(component_id=component4.id).get_state(), ControlStates.UNKNOWN)