import uuid
from datetime import datetime
from typing import List, Type

from nolabs.domain.models.common import Experiment, ExperimentId, ExperimentName
from nolabs.workflow.component import Component, JobValidationError, TInput, TOutput


class WorkflowTestsMixin:
    def seed_empty_component(self, t_input, t_output) -> Type[Component]:
        class C(Component[t_input, t_output]):
            name = 'C'

            async def setup_jobs(self):
                pass

            async def prevalidate_jobs(self) -> List[JobValidationError]:
                return []

            @property
            def _input_parameter_type(self) -> Type[TInput]:
                return t_input

            @property
            def _output_parameter_type(self) -> Type[TOutput]:
                return t_output

            async def execute(self):
                pass

        return C

    def seed_experiment(self) -> Experiment:
        e = Experiment(id=ExperimentId(uuid.uuid4()), name=ExperimentName('Test'), created_at=datetime.utcnow())
        e.save()
        return e