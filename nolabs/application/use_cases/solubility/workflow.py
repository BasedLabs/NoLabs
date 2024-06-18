import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.application.use_cases.solubility.api_models import SetupJobRequest
from nolabs.application.use_cases.solubility.use_cases import SetupJobFeature, RunJobFeature, GetJobFeature
from nolabs.domain.models.common import Protein, Experiment
from nolabs.domain.models.solubility import SolubilityJob
from nolabs.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import Component, JobValidationError


class SolubilityComponentInput(BaseModel):
    proteins: List[uuid.UUID]


class SolubilityComponentOutput(BaseModel):
    proteins_with_solubility: List[uuid.UUID]


class SolubilityComponent(Component[SolubilityComponentInput, SolubilityComponentOutput]):
    name = 'Solubility'
    description = 'Protein solubility prediction'

    async def execute(self):
        run_job_feature = RunJobFeature(api=InfrastructureDependencies.solubility_microservice())
        get_job_feature = GetJobFeature()

        for job in self.jobs:
            await run_job_feature.handle(job_id=job.id)

        items: List[uuid.UUID] = []

        for job in self.jobs:
            get_result = await get_job_feature.handle(job_id=job.id)

            for protein_id in get_result.protein_ids:
                items.append(protein_id)

        self.output = SolubilityComponentOutput(
            proteins_with_solubility=items
        )

    async def setup_jobs(self):
        setup_job_feature = SetupJobFeature()

        self.jobs = []

        for protein_id in self.input.proteins:
            protein = Protein.objects.with_id(protein_id)

            result = await setup_job_feature.handle(request=SetupJobRequest(
                experiment_id=self.experiment.id,
                protein_ids=[protein.id],
                job_name=f'Solubility {protein.name.fasta_name}'
            ))

            self.jobs.append(SolubilityJob.objects.with_id(result.job_id))

    async def jobs_setup_errors(self) -> List[JobValidationError]:
        validation_errors = []

        job: SolubilityJob
        for job in self.jobs:
            if not job.proteins:
                validation_errors.append(
                    JobValidationError(
                        job_id=job.id,
                        msg=f'No proteins'
                    )
                )

        return validation_errors

    @property
    def _input_parameter_type(self) -> Type[SolubilityComponentInput]:
        return SolubilityComponentInput

    @property
    def _output_parameter_type(self) -> Type[SolubilityComponentOutput]:
        return SolubilityComponentOutput
