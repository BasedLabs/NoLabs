import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.refined.application.use_cases.folding.api_models import FoldingBackendEnum, SetupJobRequest
from nolabs.refined.application.use_cases.folding.use_cases import SetupJobFeature, RunJobFeature, GetJobFeature
from nolabs.refined.domain.models.common import Protein, Experiment
from nolabs.refined.domain.models.folding import FoldingJob
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import PythonComponent, JobValidationError


class FoldingComponentInput(BaseModel):
    protein_ids: List[uuid.UUID]


class FoldingComponentOutput(BaseModel):
    protein_ids: List[uuid.UUID]


class FoldingComponent(PythonComponent[FoldingComponentInput, FoldingComponentOutput]):
    name = 'Folding'

    async def execute(self):
        run_job_feature = RunJobFeature(
            esmfold=InfrastructureDependencies.esmfold_microservice()
        )
        get_job_feature = GetJobFeature()

        for job in self.jobs:
            await run_job_feature.handle(job_id=job.id)

        items: List[uuid.UUID] = []

        for job in self.jobs:
            get_result = await get_job_feature.handle(job_id=job.id)

            for protein_id in get_result.proteins:
                items.append(protein_id)

        self.output = FoldingComponentOutput(
            protein_ids=items
        )

    async def setup_jobs(self):
        setup_job_feature = SetupJobFeature()

        self.jobs = []

        for protein_id in self.input.protein_ids:
            protein = Protein.objects.with_id(protein_id)

            result = await setup_job_feature.handle(request=SetupJobRequest(
                experiment_id=self.experiment.id,
                backend=FoldingBackendEnum.esmfold,
                proteins=[protein.id],
                job_id=None,
                job_name=f'Folding {protein.name.fasta_name}'
            ))

            self.jobs.append(FoldingJob.objects.with_id(result.job_id))

        self.save(cascade=True)

    async def prevalidate_jobs(self) -> List[JobValidationError]:
        validation_errors = []

        job: FoldingJob
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
    def _input_parameter_type(self) -> Type[FoldingComponentInput]:
        return FoldingComponentInput

    @property
    def _output_parameter_type(self) -> Type[FoldingComponentOutput]:
        return FoldingComponentOutput
