import uuid
from typing import List, Type

from pydantic import BaseModel

from domain.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.msa_generation.api_models import SetupJobRequest
from nolabs.application.use_cases.msa_generation.use_cases import SetupJobFeature, RunJobFeature
from nolabs.domain.models.msa import MsaGenerationJob
from nolabs.infrastructure.di import InfrastructureDependencies
from nolabs.application.workflow.component import Component, JobValidationError


class MsaGenerationInput(BaseModel):
    proteins: List[uuid.UUID]


class MsaGenerationOutput(BaseModel):
    proteins_with_msa: List[uuid.UUID]


class MsaGenerationComponent(Component[MsaGenerationInput, MsaGenerationOutput]):
    name = 'Msa generation'

    async def execute(self):
        if await self.jobs_setup_errors():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        api = InfrastructureDependencies.msa_light_microservice()
        settings = InfrastructureDependencies.msa_light_settings()
        run_job_feature = RunJobFeature(api=api,
                                        settings=settings)

        for job in self.jobs:
            await run_job_feature.handle(job_id=job.id)

        self.output = MsaGenerationOutput(
            proteins_with_msa=self.input.proteins
        )

    async def setup_jobs(self):
        self.jobs = []

        setup_job_feature = SetupJobFeature()

        for protein_id in self.input.proteins:
            job = await setup_job_feature.handle(
                request=SetupJobRequest(
                    experiment_id=self.experiment.id,
                    protein_id=protein_id
                )
            )

            job = MsaGenerationJob.objects.with_id(job.job_id)

            self.jobs.append(
                job
            )

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
    def _input_parameter_type(self) -> Type[MsaGenerationInput]:
        return MsaGenerationInput

    @property
    def _output_parameter_type(self) -> Type[MsaGenerationOutput]:
        return MsaGenerationOutput
