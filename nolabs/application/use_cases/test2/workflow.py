import uuid
from typing import Type, Optional, List, Any, ClassVar

from airflow.utils.context import Context
from pydantic import BaseModel

from nolabs.application.workflow.component import Component, TOutput, TInput
from nolabs.application.workflow.operators import SetupOperator, ExecuteJobOperator, OutputOperator
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


class Test2SetupOperator(SetupOperator):

    async def execute_async(self, context: Context) -> List[str]:
        component: Optional[Test2Component] = Component.get(self.component_id)
        experiment: Experiment = Experiment.objects.with_id(component.experiment_id)

        job_ids = []

        for i in range(component.input_value.x):
            job = Test2Job(id=JobId(uuid.uuid4()), name=JobName('test'), experiment=experiment)
            await job.save()

            job_ids.append(job.iid)

        return self.serialize_job_ids(job_ids)


class Test2ExecuteJobOperator(ExecuteJobOperator):

    async def execute_async(self, context: Context) -> Any:
        self.log.info(f'Hello there, job id is {self.job_id}')


class Test2OutputOperator(OutputOperator):
    async def execute_async(self, context: Context) -> Any:
        self.setup_output(Test2Output(x=15, y=35))


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
    def setup_operator_type(self) -> Type['SetupOperator']:
        return Test2SetupOperator

    @property
    def job_operator_type(self) -> Optional[Type['ExecuteJobOperator']]:
        return Test2ExecuteJobOperator

    @property
    def output_operator_type(self) -> Type['OutputOperator']:
        return Test2OutputOperator
