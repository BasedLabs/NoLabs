import uuid
from typing import Type, Optional, List, Any, ClassVar

from pydantic import BaseModel

from nolabs.application.workflow.component import Component, TOutput, TInput
from nolabs.application.workflow.prefect.tasks import SetupTask, ExecuteJobTask, OutputTask
from nolabs.domain.models.common import Job, JobInputError, JobId, JobName, Experiment


class Test2Input(BaseModel):
    x: int
    y: int = 10


class Test2Output(BaseModel):
    x: int
    y: int


class Test2Job(Job):

    def result_valid(self) -> bool:
        return True

    def _input_errors(self) -> List[JobInputError]:
        return []


class Test2SetupTask(SetupTask):

    async def _execute(self) -> List[uuid.UUID]:
        experiment_id = self.extra['experiment_id']

        component: Optional[Test2Component] = Component.get(self.component_id)
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        job_ids = []

        for i in range(component.input_value.x):
            job = Test2Job(id=JobId(uuid.uuid4()), name=JobName('test'), experiment=experiment)
            await job.save()

            job_ids.append(job.iid)

        return [j.value for j in job_ids]


class Test2ExecuteJobTask(ExecuteJobTask):

    async def _execute(self, job_id: uuid.UUID) -> Optional[BaseModel]:
        self.logger.info(f'Hello there, job id is {job_id}')
        return None


class Test2OutputTask(OutputTask):
    async def _execute(self) -> Optional[BaseModel]:
        self.setup_output(Test2Output(x=15, y=35))
        return None


class Test2Component(Component[Test2Input, Test2Output]):
    name: ClassVar[str] = 'test2'
    description: ClassVar[str] = 'test2 desc'

    @property
    def input_parameter_type(cls) -> Type[TInput]:
        return Test2Input

    @property
    def output_parameter_type(cls) -> Type[TOutput]:
        return Test2Output

    @property
    def setup_task_type(self) -> Type['SetupTask']:
        return Test2SetupTask

    @property
    def job_task_type(self) -> Optional[Type['ExecuteJobTask']]:
        return Test2ExecuteJobTask

    @property
    def output_task_type(self) -> Type['OutputTask']:
        return Test2OutputTask
