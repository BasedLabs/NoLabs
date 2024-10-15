import asyncio
import uuid
from typing import Type, List

import pytest
from pydantic import BaseModel

from integration.mixins import SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin
from nolabs.domain.models.common import Job, JobId, JobName
from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler
from nolabs.workflow.core.graph import GraphExecutionNode
from nolabs.workflow.core.states import ControlStates
from nolabs.infrastructure.celery_app_factory import get_celery_app
from nolabs.infrastructure.settings import settings


celery = get_celery_app()
@celery.task(name="long_running_job_test_success", bind=True, queue=settings.workflow_queue)
def long_running_job_success(bind, job_id: uuid.UUID):
    async def _():
        job: Job = Job.objects.with_id(job_id)
        job.set_name(JobName("long_running_job_test_success"))
        await job.save()
    asyncio.run(_())


class TestJobs(SeedComponentsMixin, SeedExperimentMixin, GraphTestMixin):
    @pytest.mark.asyncio
    async def test_should_successfully_run_component_with_jobs(self, prefork_celery_worker):
        j1_id = uuid.uuid4()
        j2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_component_task(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(id=JobId(j1_id), name=JobName("hello 1"), component=self.component_id)
                job2 = Job.create(id=JobId(j2_id), name=JobName("hello 2"), component=self.component_id)
                await job1.save()
                await job2.save()

                return [job1.id, job2.id]

            async def on_job_task(self, job_id: uuid.UUID):
                j: Job = Job.objects.with_id(job_id)
                j.name = JobName("Changed")
                await j.save()

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

        # act
        await graph.schedule(components=[component])
        await graph.start()
        await self.sync_until_terminal(graph=graph)

        # assert
        job1: Job = Job.objects.with_id(j1_id)
        job2: Job = Job.objects.with_id(j2_id)
        assert job1.name.value == "Changed"
        assert job2.name.value == "Changed"
        assert await graph.get_component_node(component_id=component.id).get_state() == ControlStates.SUCCESS

    @pytest.mark.asyncio
    async def test_should_fail_component_if_all_jobs_failed(self, prefork_celery_worker):
        j1_id = uuid.uuid4()
        j2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_component_task(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(id=JobId(j1_id), name=JobName("hello 1"), component=self.component_id)
                job2 = Job.create(id=JobId(j2_id), name=JobName("hello 2"), component=self.component_id)
                await job1.save()
                await job2.save()

                return [job1.id, job2.id]

            async def on_job_task(self, job_id: uuid.UUID):
                j: Job = Job.objects.with_id(job_id)
                j.name = JobName("Changed")
                await j.save()
                raise ValueError("Hello there")

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

        # act
        await graph.schedule(components=[component])
        await graph.start()
        await self.sync_until_terminal(graph=graph)

        # assert
        job1: Job = Job.objects.with_id(j1_id)
        job2: Job = Job.objects.with_id(j2_id)
        assert job1.name.value == "Changed"
        assert job2.name.value == "Changed"
        assert await graph.get_component_node(component_id=component.id).get_state() == ControlStates.FAILURE

    @pytest.mark.asyncio(loop_scope="module")
    async def test_should_successfully_run_long_running_job(self, prefork_celery_worker):
        j1_id = uuid.uuid4()
        j2_id = uuid.uuid4()

        class IO(BaseModel):
            a: int = 10

        class FlowHandler(ComponentFlowHandler):
            async def on_component_task(self, inp: IO) -> List[uuid.UUID]:
                job1 = Job.create(id=JobId(j1_id), name=JobName("hello 1"), component=self.component_id)
                job2 = Job.create(id=JobId(j2_id), name=JobName("hello 2"), component=self.component_id)
                await job1.save()
                await job2.save()

                return [job1.id, job2.id]

            async def on_job_task(self, job_id: uuid.UUID):
                j: Job = Job.objects.with_id(job_id)
                j.name = JobName("Changed")
                await j.save()
                await self.schedule_long_running(job_id=job_id, celery_task_name="long_running_job_success", input={"job_id": job_id})

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

        # act
        await graph.schedule(components=[component])
        await graph.start()
        await self.sync_until_terminal(graph=graph, timeout=10000)

        # assert
        job1: Job = Job.objects.with_id(j1_id)
        job2: Job = Job.objects.with_id(j2_id)
        assert await graph.get_component_node(component_id=component.id).get_state() == ControlStates.SUCCESS
        assert job1.name.value == "long_running_job_test_success"
        assert job2.name.value == "long_running_job_test_success"
