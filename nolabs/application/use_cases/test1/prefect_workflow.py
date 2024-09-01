import uuid
from typing import Type, Optional, List, Any, ClassVar

from airflow.utils.context import Context
from pydantic import BaseModel

from nolabs.application.workflow.component import Component, TOutput, TInput
from nolabs.application.workflow.prefect.tasks import SetupTask, ExecuteJobTask, OutputTask
from nolabs.domain.models.common import Job, JobInputError, JobId, JobName, Experiment


class Test1Input(BaseModel):
    x: int
    y: int = 10


class Test1Output(BaseModel):
    x: int
    y: int


class Test1Job(Job):

    def result_valid(self) -> bool:
        return True

    def _input_errors(self) -> List[JobInputError]:
        return []


class Test1SetupTask(SetupTask):

    async def _execute(self) -> List[uuid.UUID]:
        experiment_id = self.extra['experiment_id']

        component: Optional[Test1Component] = Component.get(self.component_id)
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        job_ids = []

        for i in range(component.input_value.x):
            job = Test1Job(id=JobId(uuid.uuid4()), name=JobName('test'), experiment=experiment)
            await job.save()

            job_ids.append(job.iid)

        return [j.value for j in job_ids]


class Test1OutputTask(OutputTask):
    async def _execute(self) -> Optional[BaseModel]:
        self.setup_output(Test1Output(x=15, y=35))
        return None


class Test1Component(Component[Test1Input, Test1Output]):
    name: ClassVar[str] = 'test1'
    description: ClassVar[str] = 'test1 desc'

    @property
    def input_parameter_type(cls) -> Type[TInput]:
        return Test1Input

    @property
    def output_parameter_type(cls) -> Type[TOutput]:
        return Test1Output

    @property
    def setup_task_type(self) -> Type['SetupTask']:
        return Test1SetupTask

    @property
    def job_task_type(self) -> Optional[Type['ExecuteJobTask']]:
        return None

    @property
    def output_task_type(self) -> Type['OutputTask']:
        return Test1OutputTask
