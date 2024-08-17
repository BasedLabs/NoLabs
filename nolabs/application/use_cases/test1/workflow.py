import uuid
from typing import Type, Optional, List, Any, ClassVar

from airflow.utils.context import Context
from pydantic import BaseModel

from nolabs.application.workflow.component import Component, TOutput, TInput
from nolabs.application.workflow.operators import SetupOperator, JobOperator, OutputOperator
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


class Test1SetupOperator(SetupOperator):

    def execute(self, context: Context) -> List[JobId]:
        component: Optional[Test1Component] = self.repository.fetch_component(self.component_id)
        experiment: Experiment = Experiment.objects.with_id(component.experiment_id)

        job_ids = []

        for i in range(component.input_value.x):
            job = Test1Job(id=JobId(uuid.uuid4()), name=JobName('test'), experiment=experiment)
            job.save()

            job_ids.append(job.iid)

        return job_ids


class Test1JobOperator(JobOperator):

    def execute(self, context: Context) -> Any:
        self.log.info(f'Hello there, job id is {self.job_id}')


class Test1OutputOperator(OutputOperator):
    def execute(self, context: Context) -> Any:
        component = self.repository.fetch_component(self.component_id)
        component.output_value = Test1Output(x=15, y=35)


class Test1Component(Component[Test1Input, Test1Output]):
    name: ClassVar[str] = 'test1'
    description: ClassVar[str] = 'test1 desc'

    @classmethod
    def input_parameter_type(cls) -> Type[TInput]:
        return Test1Input

    @classmethod
    def output_parameter_type(cls) -> Type[TOutput]:
        return Test1Output

    @property
    def setup_operator_type(self) -> Type['SetupOperator']:
        return Test1SetupOperator

    @property
    def job_operator_type(self) -> Optional[Type['JobOperator']]:
        return Test1JobOperator

    @property
    def output_operator_type(self) -> Type['OutputOperator']:
        return Test1OutputOperator
