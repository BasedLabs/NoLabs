import uuid
from typing import List, Type

from domain.exceptions import ErrorCodes, NoLabsException
from pydantic import BaseModel

from nolabs.application.use_cases.protein_design.use_cases import RunJobFeature
from nolabs.domain.models.common import JobId, JobName, Protein
from nolabs.domain.models.protein_design import ProteinDesignJob
from nolabs.infrastructure.di import InfrastructureDependencies
from application.workflow.component import Component, JobValidationError


class ProteinDesignInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class ProteinDesignOutput(BaseModel):
    binder_proteins: List[uuid.UUID]


class ProteinDesignComponent(Component[ProteinDesignInput, ProteinDesignOutput]):
    name = "Protein binder design"
    description = "Protein binder prediction using Rfdiffusion"

    async def execute(self):
        if await self.jobs_setup_errors():
            raise NoLabsException(ErrorCodes.invalid_job_input, "Jobs are not valid")

        run_job_feature = RunJobFeature(
            api=InfrastructureDependencies.protein_design_microservice()
        )

        protein_ids = []

        for job in self.jobs:
            result = await run_job_feature.handle(job_id=job.id)
            for protein_id in result.binder_ids:
                protein_ids.append(protein_id)

        self.output = ProteinDesignOutput(binder_proteins=protein_ids)

    async def setup_jobs(self):
        self.jobs = []

        for protein_id in self.input.proteins_with_pdb:
            protein = Protein.objects.with_id(protein_id)

            if not protein.pdb_content:
                raise NoLabsException(
                    ErrorCodes.protein_pdb_is_empty, "Protein pdb content is undefined"
                )

            job_id = JobId(uuid.uuid4())
            job_name = JobName(f"Protein binder design for protein {protein.name}")

            job = ProteinDesignJob(id=job_id, name=job_name, experiment=self.experiment)

            job.set_protein(protein=protein)

            job.save()

            self.jobs.append(job)

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
    def _input_parameter_type(self) -> Type[ProteinDesignInput]:
        return ProteinDesignInput

    @property
    def _output_parameter_type(self) -> Type[ProteinDesignOutput]:
        return ProteinDesignOutput
