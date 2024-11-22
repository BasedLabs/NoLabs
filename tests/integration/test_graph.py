import asyncio
import uuid
from multiprocessing import Value
from typing import List, Optional, Type

from asgiref.sync import async_to_sync
from pydantic import BaseModel

from nolabs.domain.models.common import Job, JobId, JobName
from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.infrastructure.settings import settings
from nolabs.workflow.core.component import Component
from nolabs.workflow.core.flow import ComponentFlowHandler
from nolabs.workflow.core.graph import Graph
from nolabs.workflow.core.states import ControlStates
from tests.integration.mixins import (
    GraphTestMixin,
    SeedComponentsMixin,
    SeedExperimentMixin,
)
from tests.integration.setup import GlobalSetup


class TestGraph(GlobalSetup, SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin):
    async def test_should_complete_graph_if_component_failed(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
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
            await graph.get_component_node(component.id).get_state(),
            ControlStates.FAILURE,
        )

    async def test_should_run_linear_workflow(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_finish(
                self, inp: IO, job_ids: List[uuid.UUID]
            ) -> Optional[IO]:
                return IO(a=15)

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
        component1 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component2 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        self.seed_mappings(component2, previous_components=[(component1, ["a"], ["a"])])
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component1, component2])
        await graph.schedule(component_ids=[component1.id, component2.id])
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component1.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component2.id).get_state(),
            ControlStates.SUCCESS,
        )

    async def test_should_run_diamond_workflow(self):
        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_finish(
                self, inp: IO, job_ids: List[uuid.UUID]
            ) -> Optional[IO]:
                return IO(a=15)

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
        component1 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component2 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component3 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component4 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        self.seed_mappings(component2, previous_components=[(component1, ["a"], ["a"])])
        self.seed_mappings(component3, previous_components=[(component1, ["a"], ["a"])])
        self.seed_mappings(
            component4,
            previous_components=[
                (component2, ["a"], ["a"]),
                (component3, ["a"], ["a"]),
            ],
        )
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(
            components=[component1, component2, component3, component4]
        )
        await graph.schedule(
            component_ids=[component1.id, component2.id, component3.id, component4.id]
        )
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component1.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component2.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component3.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component4.id).get_state(),
            ControlStates.SUCCESS,
        )

    async def test_should_not_propagate_execution_if_component_failed(self):
        component2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_finish(
                self, inp: IO, job_ids: List[uuid.UUID]
            ) -> Optional[IO]:
                if self.component_id == component2_id:
                    raise ValueError("Component failed to propagate")
                return IO(a=15)

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
        component1 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component2 = self.seed_component(
            experiment_id=experiment_id,
            component_type=MockComponent,
            component_id=component2_id,
        )
        component3 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component4 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        self.seed_mappings(component2, previous_components=[(component1, ["a"], ["a"])])
        self.seed_mappings(component3, previous_components=[(component1, ["a"], ["a"])])
        self.seed_mappings(
            component4,
            previous_components=[
                (component2, ["a"], ["a"]),
                (component3, ["a"], ["a"]),
            ],
        )
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(
            components=[component1, component2, component3, component4]
        )
        await graph.schedule(
            component_ids=[component1.id, component2.id, component3.id, component4.id]
        )
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component1.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component2.id).get_state(),
            ControlStates.FAILURE,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component3.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component4.id).get_state(),
            ControlStates.CANCELLED,
        )

    async def test_can_schedule_separate_components(self):
        component2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_finish(
                self, inp: IO, job_ids: List[uuid.UUID]
            ) -> Optional[IO]:
                return IO(a=15)

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
        component1 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component2 = self.seed_component(
            experiment_id=experiment_id,
            component_type=MockComponent,
            component_id=component2_id,
        )
        self.seed_mappings(component2, previous_components=[(component1, ["a"], ["a"])])
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component1, component2])
        await graph.schedule(component_ids=[component1.id])
        await graph.sync()

        await graph.schedule(component_ids=[component2.id])
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component1.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component2.id).get_state(),
            ControlStates.SUCCESS,
        )

    async def test_component_is_not_called_second_time(self):
        component2_id = uuid.uuid4()
        shared_counter = Value("i", 0)

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_finish(
                self, inp: IO, job_ids: List[uuid.UUID]
            ) -> Optional[IO]:
                shared_counter.value += 1
                return IO(a=15)

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
        component1 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component2 = self.seed_component(
            experiment_id=experiment_id,
            component_type=MockComponent,
            component_id=component2_id,
        )
        self.seed_mappings(component2, previous_components=[(component1, ["a"], ["a"])])
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component1, component2])
        await graph.schedule(component_ids=[component1.id])
        await graph.sync()

        await graph.schedule(component_ids=[component2.id])
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component1.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component2.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(shared_counter.value, 2)

    async def test_can_schedule_another_component_while_first_running(self):
        task_name = str(uuid.uuid4())
        shared_counter = Value("i", 0)

        class IO(BaseModel):
            a: int = 10

        celery = get_celery_app()

        def task(bind, job_id: uuid.UUID):
            async def _():
                print("Cycle")
                shared_counter.value += 1
                while shared_counter.value != 2:
                    await asyncio.sleep(0.1)
                print("Exit")

            async_to_sync(_)()

        celery.task(task, name=task_name, bind=True, queue=settings.workflow_queue)

        class FlowHandler(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(
                    id=JobId(uuid.uuid4()),
                    name=JobName("hello 1"),
                    component=self.component_id,
                )
                await job1.save()

                print("Job created")

                return [job1.id]

            async def on_job_start(self, job_id: uuid.UUID):
                j: Job = Job.objects.with_id(job_id)
                j.name = JobName("Changed")
                await j.save()
                await self.schedule(
                    job_id=job_id,
                    celery_queue=settings.workflow_queue,
                    celery_task_name=task_name,
                    input={"job_id": job_id}
                )

            async def on_finish(
                self, inp: IO, job_ids: List[uuid.UUID]
            ) -> Optional[IO]:
                return IO(a=145)

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
        component1 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        component2 = self.seed_component(
            experiment_id=experiment_id, component_type=MockComponent
        )
        graph = Graph(experiment_id=experiment_id)

        self.spin_up_celery()

        # act
        await graph.set_components_graph(components=[component1, component2])
        await graph.schedule(component_ids=[component1.id])
        await graph.sync()

        await graph.schedule(component_ids=[component2.id])

        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            await graph.get_component_node(component_id=component1.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(
            await graph.get_component_node(component_id=component2.id).get_state(),
            ControlStates.SUCCESS,
        )
