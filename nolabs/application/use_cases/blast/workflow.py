import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.application.use_cases.blast.api_models import SetupJobRequest
from nolabs.application.use_cases.blast.use_cases import SetupJobFeature, RunJobFeature, GetJobFeature
from nolabs.domain.models.common import Protein
from nolabs.domain.models.blast import BlastJob
from nolabs.infrastructure.di import InfrastructureDependencies
from nolabs.application.workflow.component import Component, JobValidationError


class BlastComponentInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]
    descriptions: int = 10
    alignments: int = 10
    hitlist_size: int = 10
    expect: int = 10



class BlastComponentOutput(BaseModel):
    hits: List[str]


class BlastComponent(Component[BlastComponentInput, BlastComponentOutput]):
    name = 'Blast'
    description = 'Finds similar sequences for proteins and nucleotides.'

    async def execute(self):
        run_job_feature = RunJobFeature(blast=InfrastructureDependencies.blast_query_microservice())
        get_job_feature = GetJobFeature()

        for job in self.jobs:
            await run_job_feature.handle(job_id=job.id)

        items: List[str] = []

        for job in self.jobs:
            get_result = await get_job_feature.handle(job_id=job.id)

            for result in get_result.result:
                for hit in result.hits:
                    items.append(
                        hit.id
                    )

        self.output = BlastComponentOutput(
            hits=items
        )

    async def setup_jobs(self):
        setup_job_feature = SetupJobFeature()

        self.jobs = []

        for protein_id in self.input.proteins_with_fasta:
            protein: Protein = Protein.objects.with_id(protein_id)

            result = await setup_job_feature.handle(request=SetupJobRequest(
                experiment_id=self.experiment.id,
                protein_id=protein.iid.value,
                descriptions=self.input.descriptions,
                alignments=self.input.alignments,
                hitlist_size=self.input.hitlist_size,
                expect=self.input.expect,
                job_id=None,
                job_name=f'Blast {protein.name.fasta_name}'
            ))

            self.jobs.append(BlastJob.objects.with_id(result.job_id))

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
    def _input_parameter_type(self) -> Type[BlastComponentInput]:
        return BlastComponentInput

    @property
    def _output_parameter_type(self) -> Type[BlastComponentOutput]:
        return BlastComponentOutput
