import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.conformations.use_cases import RunJobFeature
from nolabs.domain.models.common import Protein, JobId, JobName
from nolabs.domain.models.conformations import ConformationsJob
from nolabs.infrastructure.di import InfrastructureDependencies
from nolabs.application.workflow.component import Component, JobValidationError


class ConformationInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class ConformationOutput(BaseModel):
    proteins_with_conformations: List[uuid.UUID]


class ConformationComponent(Component[ConformationInput, ConformationOutput]):
    name = 'Conformations'
    description = 'Protein molecular dynamics'

    async def execute(self):
        if await self.jobs_setup_errors():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        run_job_feature = RunJobFeature(api=InfrastructureDependencies.conformations_microservice())

        protein_ids = []

        for job in self.jobs:
            result = await run_job_feature.handle(job_id=job.id)
            protein_ids.append(result.protein_id)

        self.output = ConformationOutput(proteins_with_conformations=protein_ids)

    async def setup_jobs(self):
        self.jobs = []

        for protein_id in self.input.proteins_with_pdb:
            protein = Protein.objects.with_id(protein_id)

            job_id = JobId(uuid.uuid4())
            job_name = JobName(f'Conformations for protein {protein.name}')

            job = ConformationsJob(
                id=job_id,
                name=job_name,
                experiment=self.experiment
            )

            job.set_protein(
                protein=protein
            )

            await job.save()

            self.jobs.append(job)

    async def jobs_setup_errors(self) -> List[JobValidationError]:
        jobs_errors = []

        for job in self.jobs:
            errors = job.input_errors()
            if errors:
                jobs_errors.append(
                    JobValidationError(
                        job_id=job.id,
                        msg=', '.join([err.message for err in errors])
                    )
                )

        return jobs_errors

    @property
    def _input_parameter_type(self) -> Type[ConformationInput]:
        return ConformationInput

    @property
    def _output_parameter_type(self) -> Type[ConformationOutput]:
        return ConformationOutput
