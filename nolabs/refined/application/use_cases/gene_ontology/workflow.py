import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.refined.application.use_cases.gene_ontology.api_models import SetupJobRequest
from nolabs.refined.application.use_cases.gene_ontology.use_cases import SetupJobFeature, RunJobFeature, GetJobFeature
from nolabs.refined.domain.models.common import Protein, Experiment
from nolabs.refined.domain.models.gene_ontology import GeneOntologyJob
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import Component, JobValidationError


class GeneOntologyComponentInput(BaseModel):
    proteins: List[uuid.UUID]


class GeneOntologyComponentOutput(BaseModel):
    proteins_with_gene_ontology: List[uuid.UUID]


class GeneOntologyComponent(Component[GeneOntologyComponentInput, GeneOntologyComponentOutput]):
    name = 'Gene ontology'
    description = 'Protein gene ontology prediction'

    async def execute(self):
        run_job_feature = RunJobFeature(api=InfrastructureDependencies.gene_ontology_microservice())
        get_job_feature = GetJobFeature()

        for job in self.jobs:
            await run_job_feature.handle(job_id=job.id)

        items: List[uuid.UUID] = []

        for job in self.jobs:
            get_result = await get_job_feature.handle(job_id=job.id)

            for protein_id in get_result.protein_ids:
                items.append(protein_id)

        self.output = GeneOntologyComponentOutput(
            proteins_with_gene_ontology=items
        )

    async def setup_jobs(self):
        setup_job_feature = SetupJobFeature()

        self.jobs = []

        for protein_id in self.input.proteins:
            protein = Protein.objects.with_id(protein_id)

            result = await setup_job_feature.handle(request=SetupJobRequest(
                experiment_id=self.experiment.id,
                proteins=[protein.id],
                job_name=f'GeneOntology {protein.name.fasta_name}'
            ))

            self.jobs.append(GeneOntologyJob.objects.with_id(result.job_id))

    async def prevalidate_jobs(self) -> List[JobValidationError]:
        validation_errors = []

        job: GeneOntologyJob
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
    def _input_parameter_type(self) -> Type[GeneOntologyComponentInput]:
        return GeneOntologyComponentInput

    @property
    def _output_parameter_type(self) -> Type[GeneOntologyComponentOutput]:
        return GeneOntologyComponentOutput
