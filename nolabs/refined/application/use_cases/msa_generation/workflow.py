import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.msa_generation.api_models import SetupJobRequest
from nolabs.refined.application.use_cases.msa_generation.use_cases import SetupJobFeature, RunJobFeature
from nolabs.refined.domain.models.msa import MsaGenerationJob
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import Component, JobValidationError


class MsaGenerationInput(BaseModel):
    proteins: List[uuid.UUID]


class MsaGenerationOutput(BaseModel):
    proteins_with_msa: List[uuid.UUID]


class MsaGenerationComponent(Component[MsaGenerationInput, MsaGenerationOutput]):
    name = 'Msa generation'

    async def execute(self):
        if not self.prevalidate_jobs():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        run_job_feature = RunJobFeature(api=InfrastructureDependencies.protein_design_microservice())

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

    async def prevalidate_jobs(self) -> List[JobValidationError]:
        validation_errors = []

        job: MsaGenerationJob
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
    def _input_parameter_type(self) -> Type[MsaGenerationInput]:
        return MsaGenerationInput

    @property
    def _output_parameter_type(self) -> Type[MsaGenerationOutput]:
        return MsaGenerationOutput
