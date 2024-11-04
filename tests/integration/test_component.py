import asyncio
import uuid
from typing import List, Optional, Type

from pydantic import BaseModel

from nolabs.domain.models.common import Job, JobId, JobName
from nolabs.workflow.core.component import Component, TInput, TOutput
from nolabs.workflow.core.flow import ComponentFlowHandler
from nolabs.workflow.core.graph import Graph
from nolabs.workflow.core.states import ControlStates
from tests.integration.mixins import (
    GraphTestMixin,
    SeedComponentsMixin,
    SeedExperimentMixin,
)
from tests.integration.setup import GlobalSetup


class TestComponent(
    GlobalSetup, SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin
):
    async def test_should_successfully_run_component(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler): ...

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[IO]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return FlowHandler

            @property
            def output_parameter_type(self) -> Type[IO]:
                return IO

        # arrange
        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component])
        await graph.schedule(component_ids=[component.id])
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_state(),
            ControlStates.SUCCESS,
        )

    async def test_should_fail_component_on_main_task_failure(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_component_task(self, inp: TInput) -> List[uuid.UUID]:
                raise ValueError("Hello")

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[IO]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return FlowHandler

            @property
            def output_parameter_type(self) -> Type[IO]:
                return IO

        # arrange
        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component])
        await graph.schedule(component_ids=[component.id])
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_state(),
            ControlStates.FAILURE,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_message(),
            "Hello",
        )

    async def test_should_fail_component_on_complete_task_failure(self):
        # arrange

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_completion(
                self, inp: IO, job_ids: List[uuid.UUID]
            ) -> Optional[TOutput]:
                raise ValueError("Hello")

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[IO]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return FlowHandler

            @property
            def output_parameter_type(self) -> Type[IO]:
                return IO

        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component])
        await graph.schedule(component_ids=[component.id])
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_state(),
            ControlStates.FAILURE,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_message(),
            "Hello",
        )

    async def test_should_fail_component_on_all_jobs_failure(self):
        # arrange

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_component_task(self, inp: TInput) -> List[uuid.UUID]:
                job1 = Job.create(
                    id=JobId(uuid.uuid4()),
                    name=JobName("hello 1"),
                    component=self.component_id,
                )
                job2 = Job.create(
                    id=JobId(uuid.uuid4()),
                    name=JobName("hello 2"),
                    component=self.component_id,
                )
                await job1.save()
                await job2.save()

                return [job1.id, job2.id]

            def on_job_task(self, job_id: uuid.UUID):
                raise ValueError("Exception")

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[IO]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return FlowHandler

            @property
            def output_parameter_type(self) -> Type[IO]:
                return IO

        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component])
        await graph.schedule(component_ids=[component.id])
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_state(),
            ControlStates.FAILURE,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_message(),
            "All jobs failed",
        )

    async def test_should_cancel_component(self):
        # arrange

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_component_task(self, inp: IO) -> List[uuid.UUID]:
                await asyncio.sleep(1000)
                return []

            def on_job_task(self, job_id: uuid.UUID):
                raise ValueError("Exception")

        class MockComponent(Component[IO, IO], ComponentFlowHandler):
            name = "a"

            @property
            def input_parameter_type(self) -> Type[IO]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return FlowHandler

            @property
            def output_parameter_type(self) -> Type[IO]:
                return IO

        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component])
        await graph.schedule(component_ids=[component.id])
        await graph.sync()
        await asyncio.sleep(1)
        await graph.cancel()
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_state(),
            ControlStates.CANCELLED,
        )
