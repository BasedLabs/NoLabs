import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.conformations.use_cases import RunJobFeature
from nolabs.refined.domain.models.common import Protein, JobId, JobName
from nolabs.refined.domain.models.conformations import ConformationsJob
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import Component, JobValidationError


class ConformationInput(BaseModel):
    proteins: List[uuid.UUID]


class ConformationOutput(BaseModel):
    proteins_with_conformations: List[uuid.UUID]


class ConformationComponent(Component[ConformationInput, ConformationOutput]):
    name = 'Conformations'

    async def execute(self):
        if not self.prevalidate_jobs():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        run_job_feature = RunJobFeature(api=InfrastructureDependencies.conformations_microservice())

        protein_ids = []

        for job in self.jobs:
            result = await run_job_feature.handle(job_id=job.id)
            protein_ids.append(result.protein_id)

        self.output = ConformationOutput(proteins_with_conformations=protein_ids)

    async def setup_jobs(self):
        self.jobs = []

        for protein_id in self.input.proteins:
            protein = Protein.objects.with_id(protein_id)

            job_id = JobId(uuid.uuid4())
            job_name = JobName(f'Conformations for protein {protein.name}')

            job = ConformationsJob(
                id=job_id,
                name=job_name,
                experiment=self.experiment
            )

            job.save()

            self.jobs.append(job)

    async def prevalidate_jobs(self) -> List[JobValidationError]:
        validation_errors = []

        job: ConformationsJob
        for job in self.jobs:
            if not job.input_valid():
                validation_errors.append(
                    JobValidationError(
                        job_id=job.id,
                        msg=f'Job input is invalid'
                    )
                )

        return validation_errors

    @property
    def _input_parameter_type(self) -> Type[ConformationInput]:
        return ConformationInput

    @property
    def _output_parameter_type(self) -> Type[ConformationOutput]:
        return ConformationOutput