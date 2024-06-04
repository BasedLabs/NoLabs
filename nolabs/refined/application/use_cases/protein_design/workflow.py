import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.protein_design.use_cases import RunJobFeature
from nolabs.refined.domain.models.common import Protein, JobId, JobName
from nolabs.refined.domain.models.protein_design import ProteinDesignJob
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import Component, JobValidationError


class ProteinDesignInput(BaseModel):
    proteins: List[uuid.UUID]


class ProteinDesignOutput(BaseModel):
    binder_proteins: List[uuid.UUID]


class ProteinDesignComponent(Component[ProteinDesignInput, ProteinDesignOutput]):
    name = 'Protein binder design'

    async def execute(self):
        if not self.prevalidate_jobs():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        run_job_feature = RunJobFeature(api=InfrastructureDependencies.protein_design_microservice())

        protein_ids = []

        for job in self.jobs:
            result = await run_job_feature.handle(job_id=job.id)
            for protein_id in result.binder_ids:
                protein_ids.append(protein_id)

        self.output = ProteinDesignOutput(binder_proteins=protein_ids)

    async def setup_jobs(self):
        self.jobs = []

        for protein_id in self.input.proteins:
            protein = Protein.objects.with_id(protein_id)

            job_id = JobId(uuid.uuid4())
            job_name = JobName(f'Protein binder design for protein {protein.name}')

            job = ProteinDesignJob(
                id=job_id,
                name=job_name,
                experiment=self.experiment
            )

            job.save()

            self.jobs.append(job)

    async def prevalidate_jobs(self) -> List[JobValidationError]:
        validation_errors = []

        job: ProteinDesignJob
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
    def _input_parameter_type(self) -> Type[ProteinDesignInput]:
        return ProteinDesignInput

    @property
    def _output_parameter_type(self) -> Type[ProteinDesignOutput]:
        return ProteinDesignOutput
