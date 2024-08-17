import uuid
from typing import Type, Optional, List, Any, ClassVar

from airflow.utils.context import Context
from pydantic import BaseModel

from nolabs.application.workflow.core.component import Component, TOutput, TInput
from nolabs.application.workflow.core.operators import SetupOperator, JobOperator, OutputOperator
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

    def execute(self, context: Context) -> List[JobId]:
        component: Optional[Test2Component] = self.repository.fetch_component(self.component_id)
        experiment: Experiment = Experiment.objects.with_id(component.experiment_id)

        job_ids = []

        for i in range(component.input_value.x):
            job = Test2Job(id=JobId(uuid.uuid4()), name=JobName('test'), experiment=experiment)
            job.save()

            job_ids.append(job.iid)

        return job_ids


class Test2JobOperator(JobOperator):

    def execute(self, context: Context) -> Any:
        self.log.info(f'Hello there, job id is {self.job_id}')


class Test2OutputOperator(OutputOperator):
    def execute(self, context: Context) -> Any:
        component = self.repository.fetch_component(self.component_id)
        component.output_value = Test2Output(x=15, y=35)


class Test2Component(Component[Test2Input, Test2Output]):
    name: ClassVar[str] = 'test2'
    description: ClassVar[str] = 'test2 desc'

    @classmethod
    def input_parameter_type(cls) -> Type[TInput]:
        return Test2Input

    @classmethod
    def output_parameter_type(cls) -> Type[TOutput]:
        return Test2Output

    @property
    def setup_operator_type(self) -> Type['SetupOperator']:
        return Test2SetupOperator

    @property
    def job_operator_type(self) -> Optional[Type['JobOperator']]:
        return Test2JobOperator

    @property
    def output_operator_type(self) -> Type['OutputOperator']:
        return Test2OutputOperator
