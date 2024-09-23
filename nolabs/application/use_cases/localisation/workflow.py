import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.application.use_cases.localisation.api_models import \
    SetupJobRequest
from nolabs.application.use_cases.localisation.use_cases import (
    GetJobFeature, RunJobFeature, SetupJobFeature)
from nolabs.domain.models.common import Protein
from nolabs.domain.models.localisation import LocalisationJob
from nolabs.infrastructure.di import InfrastructureDependencies
from nolabs.application.workflow.component import Component, JobValidationError


class LocalisationComponentInput(BaseModel):
    proteins: List[uuid.UUID]


class LocalisationComponentOutput(BaseModel):
    proteins: List[uuid.UUID]


class LocalisationComponent(
    Component[LocalisationComponentInput, LocalisationComponentOutput]
):
    name = "Localisation"
    description = "Protein localisation prediction"

    async def execute(self):
        run_job_feature = RunJobFeature(
            api=InfrastructureDependencies.localisation_microservice()
        )
        get_job_feature = GetJobFeature()

        for job in self.jobs:
            await run_job_feature.handle(job_id=job.id)

        items: List[uuid.UUID] = []

        for job in self.jobs:
            get_result = await get_job_feature.handle(job_id=job.id)

            for protein_id in get_result.protein_ids:
                items.append(protein_id)

        self.output = LocalisationComponentOutput(proteins=items)

    async def setup_jobs(self):
        setup_job_feature = SetupJobFeature()

        self.jobs = []

        for protein_id in self.input.proteins:
            protein = Protein.objects.with_id(protein_id)

            result = await setup_job_feature.handle(
                request=SetupJobRequest(
                    experiment_id=self.experiment.id,
                    protein_ids=[protein.id],
                    job_id=None,
                    job_name=f"Localisation {protein.name.fasta_name}",
                )
            )

            self.jobs.append(LocalisationJob.objects.with_id(result.job_id))

    async def jobs_setup_errors(self) -> List[JobValidationError]:
        jobs_errors = []

        for job in self.jobs:
            errors = job.input_errors()
            if errors:
                jobs_errors.append(
                    JobValidationError(
                        job_id=job.id, msg=", ".join([err.message for err in errors])
                    )
                )

        return jobs_errors

    @property
    def _input_parameter_type(self) -> Type[LocalisationComponentInput]:
        return LocalisationComponentInput

    @property
    def _output_parameter_type(self) -> Type[LocalisationComponentOutput]:
        return LocalisationComponentOutput
