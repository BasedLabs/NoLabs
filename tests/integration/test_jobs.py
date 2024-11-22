import uuid
from multiprocessing import Value
from typing import Any, Dict, List, Optional, Type

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


class TestJobs(GlobalSetup, SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin):
    async def test_should_complete_without_long_running_job(self):
        j1_id = uuid.uuid4()
        j2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(
                    id=JobId(j1_id),
                    name=JobName("hello 1"),
                    component=self.component_id,
                )
                job2 = Job.create(
                    id=JobId(j2_id),
                    name=JobName("hello 2"),
                    component=self.component_id,
                )
                await job1.save()
                await job2.save()

                return [job1.id, job2.id]

            async def on_job_start(self, job_id: uuid.UUID):
                j: Job = Job.objects.with_id(job_id)
                j.name = JobName("Changed")
                await j.save()

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

    async def test_should_successfully_run_long_running_job(self):
        celery = get_celery_app()

        def long_running_job_success(bind, job_id: uuid.UUID):
            async def _():
                print("Long running success")
                job: Job = Job.objects.with_id(job_id)
                job.set_name(JobName("long_running_job_test_success"))
                await job.save()

            async_to_sync(_)()

        celery.task(
            long_running_job_success,
            name="long_running_job_test_success2",
            bind=True,
            queue=settings.workflow_queue,
        )

        j1_id = uuid.uuid4()
        j2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(
                    id=JobId(j1_id),
                    name=JobName("hello 1"),
                    component=self.component_id,
                )
                job2 = Job.create(
                    id=JobId(j2_id),
                    name=JobName("hello 2"),
                    component=self.component_id,
                )
                await job1.save()
                await job2.save()

                return [job1.id, job2.id]

            async def on_job_start(self, job_id: uuid.UUID):
                j: Job = Job.objects.with_id(job_id)
                j.name = JobName("Changed")
                await j.save()
                await self.schedule(
                    job_id=job_id,
                    celery_queue=settings.workflow_queue,
                    celery_task_name="long_running_job_test_success2",
                    input={"job_id": job_id},
                )

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
        job1 = Job.objects.get(id=j1_id)
        job2 = Job.objects.get(id=j2_id)

        self.assertEqual(
            await graph.get_component_node(component_id=component.id).get_state(),
            ControlStates.SUCCESS,
        )
        self.assertEqual(job1.name.value, "long_running_job_test_success")
        self.assertEqual(job2.name.value, "long_running_job_test_success")

    async def test_should_fail_on_longrunning_fail(self):
        celery = get_celery_app()

        def long_running_job_test_failed(bind, job_id: uuid.UUID):
            async def _():
                raise ValueError("Hello")

            async_to_sync(_)()

        celery.task(
            long_running_job_test_failed,
            name="long_running_job_test_failed",
            bind=True,
            queue=settings.workflow_queue,
        )

        j1_id = uuid.uuid4()
        j2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(
                    id=JobId(j1_id),
                    name=JobName("hello 1"),
                    component=self.component_id,
                )
                job2 = Job.create(
                    id=JobId(j2_id),
                    name=JobName("hello 2"),
                    component=self.component_id,
                )
                await job1.save()
                await job2.save()

                return [job1.id, job2.id]

            async def on_job_start(self, job_id: uuid.UUID):
                j: Job = Job.objects.with_id(job_id)
                j.name = JobName("Changed")
                await j.save()
                await self.schedule(
                    job_id=job_id,
                    celery_queue=settings.workflow_queue,
                    celery_task_name="long_running_job_test_failed",
                    input={"job_id": job_id},
                )

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
            await graph.get_job_node(
                component_id=component.id, job_id=j1_id
            ).get_message(),
            "Hello",
        )
        self.assertEqual(
            await graph.get_job_node(
                component_id=component.id, job_id=j2_id
            ).get_message(),
            "Hello",
        )

    async def test_should_fail_on_completed_failure(self):
        celery = get_celery_app()

        def long_running_job_test_success(bind, job_id: uuid.UUID):
            print("ok")

        celery.task(
            long_running_job_test_success,
            name="long_running_job_test_success",
            bind=True,
            queue=settings.workflow_queue,
        )

        j1_id = uuid.uuid4()
        j2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(
                    id=JobId(j1_id),
                    name=JobName("hello 1"),
                    component=self.component_id,
                )
                job2 = Job.create(
                    id=JobId(j2_id),
                    name=JobName("hello 2"),
                    component=self.component_id,
                )
                await job1.save()
                await job2.save()

                return [job1.id, job2.id]

            async def on_job_start(self, job_id: uuid.UUID):
                j: Job = Job.objects.with_id(job_id)
                j.name = JobName("Changed")
                await j.save()
                await self.schedule(
                    job_id=job_id,
                    celery_queue=settings.workflow_queue,
                    celery_task_name="long_running_job_test_success",
                    input={"job_id": job_id},
                )

            async def on_job_finish(
                    self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
            ):
                raise ValueError("Failed completion")

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
            await graph.get_job_node(
                component_id=component.id, job_id=j1_id
            ).get_message(),
            "Failed completion",
        )
        self.assertEqual(
            await graph.get_job_node(
                component_id=component.id, job_id=j2_id
            ).get_message(),
            "Failed completion",
        )

    async def test_should_pass_data_from_long_running_job(self):
        celery = get_celery_app()

        def lrj(bind, job_id: uuid.UUID):
            async def _():
                return {"a": 10}

            return async_to_sync(_)()

        task_name = str(uuid.uuid4())
        celery.task(lrj, name=task_name, bind=True, queue=settings.workflow_queue)

        j1_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(
                    id=JobId(j1_id),
                    name=JobName("hello 1"),
                    component=self.component_id,
                )
                await job1.save()

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

            async def on_job_finish(
                    self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
            ):
                if not long_running_output or not long_running_output.get("a"):
                    raise ValueError("Output is empty")

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

    async def test_should_not_rerun_jobs_after_next_component_failed(self):
        called_counter = Value("i", 0)
        component1_id = uuid.uuid4()
        component2_id = uuid.uuid4()
        j1 = Job.create(
            id=JobId(uuid.uuid4()),
            name=JobName("hello 1"),
            component=component1_id,
        )
        j2 = Job.create(
            id=JobId(uuid.uuid4()),
            name=JobName("hello 2"),
            component=component2_id,
        )

        class IO(BaseModel):
            a: int = 10

        class Fh1(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
                return [j1.id]

            async def on_job_start(self, job_id: uuid.UUID):
                called_counter.value += 1

        class Fh2(ComponentFlowHandler):
            async def on_start(self, inp: IO) -> List[uuid.UUID]:
                return [j2.id]

            async def on_job_finish(
                    self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
            ):
                raise ValueError("Expected")

        class MockComponent1(Component[IO, IO], ComponentFlowHandler):
            name = "MockComponent1"

            @property
            def input_parameter_type(self) -> Type[IO]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return Fh1

            @property
            def output_parameter_type(self) -> Type[IO]:
                return IO

        class MockComponent2(Component[IO, IO], ComponentFlowHandler):
            name = "MockComponent2"

            @property
            def input_parameter_type(self) -> Type[IO]:
                return IO

            @property
            def component_flow_type(self) -> Type["ComponentFlowHandler"]:
                return Fh2

            @property
            def output_parameter_type(self) -> Type[IO]:
                return IO

        # arrange
        experiment_id = uuid.uuid4()
        await self.seed_experiment(id=experiment_id)
        component1 = self.seed_component(
            experiment_id=experiment_id,
            component_type=MockComponent1,
            component_id=component1_id
        )
        component2 = self.seed_component(
            experiment_id=experiment_id,
            component_type=MockComponent2,
            component_id=component2_id
        )
        self.seed_mappings(component2, previous_components=[(component1, ["a"], ["a"])])
        graph = Graph(experiment_id=experiment_id)
        self.spin_up_celery()
        await j1.save()
        await j2.save()

        # act
        await graph.set_components_graph(components=[component1, component2])
        await graph.schedule(component_ids=[component1.id, component2.id])
        await self.spin_up_sync(graph=graph)

        await graph.schedule(component_ids=[component1.id, component2.id])
        await self.spin_up_sync(graph=graph)

        # assert
        self.assertEqual(
            1,
            called_counter.value
        )
