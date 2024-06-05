import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.refined.application.use_cases.diffdock.api_models import SetupJobRequest
from nolabs.refined.application.use_cases.diffdock.use_cases import SetupJobFeature, RunJobFeature, GetJobFeature
from nolabs.refined.domain.models.common import Protein, Ligand
from nolabs.refined.domain.models.diffdock import DiffDockBindingJob
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import Component, JobValidationError


class DiffDockComponentInput(BaseModel):
    proteins: List[uuid.UUID]
    ligands: List[uuid.UUID]


class DiffDockComponentOutput(BaseModel):
    protein_complexes: List[uuid.UUID]


class DiffDockComponent(Component[DiffDockComponentInput, DiffDockComponentOutput]):
    name = 'DiffDock'

    async def execute(self):
        run_job_feature = RunJobFeature(diffdock=InfrastructureDependencies.diffdock_microservice())
        get_job_feature = GetJobFeature()

        for job in self.jobs:
            await run_job_feature.handle(job_id=job.id)

        items: List[uuid.UUID] = []

        for job in self.jobs:
            get_result = await get_job_feature.handle(job_id=job.id)

            for result in get_result.result:
                items.append(
                    result.complex_id
                )

        self.output = DiffDockComponentOutput(
            protein_complexes=items
        )

    async def setup_jobs(self):
        setup_job_feature = SetupJobFeature()

        self.jobs = []

        for protein_id in self.input.proteins:
            for ligand_id in self.input.ligands:
                protein: Protein = Protein.objects.with_id(protein_id)
                ligand: Ligand = Ligand.objects.with_id(ligand_id)

                result = await setup_job_feature.handle(request=SetupJobRequest(
                    experiment_id=self.experiment.id,
                    protein_id=protein.iid.value,
                    ligand_id=ligand.iid.value,
                    job_id=None,
                    job_name=f'DiffDock {protein.name.fasta_name}'
                ))

                self.jobs.append(DiffDockBindingJob.objects.with_id(result.job_id))

    async def prevalidate_jobs(self) -> List[JobValidationError]:
        validation_errors = []

        job: DiffDockBindingJob
        for job in self.jobs:
            if not job.input_valid():
                validation_errors.append(
                    JobValidationError(
                        job_id=job.id,
                        msg=f'Jobs inputs are invalid'
                    )
                )

        return validation_errors

    @property
    def _input_parameter_type(self) -> Type[DiffDockComponentInput]:
        return DiffDockComponentInput

    @property
    def _output_parameter_type(self) -> Type[DiffDockComponentOutput]:
        return DiffDockComponentOutput