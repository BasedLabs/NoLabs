import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.binding_pockets.api_models import SetupJobRequest
from nolabs.application.use_cases.binding_pockets.use_cases import RunJobFeature, SetupJobFeature
from nolabs.domain.models.pocket_prediction import PocketPredictionJob
from nolabs.domain.models.common import Protein, JobId, JobName
from nolabs.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import Component, JobValidationError


class BindingPocketPredictionInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class BindingPocketPredictionOutput(BaseModel):
    proteins_with_binding_pockets: List[uuid.UUID]


class BindingPocketPredictionComponent(Component[BindingPocketPredictionInput, BindingPocketPredictionOutput]):
    name = 'Binding pockets'
    description = 'Protein binding pockets prediction'

    async def execute(self):
        if not self.prevalidate_jobs():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        run_job_feature = RunJobFeature(api=InfrastructureDependencies.p2rank_microservice())

        protein_ids = []

        for job in self.jobs:
            result = await run_job_feature.handle(job_id=job.id)
            protein_ids.append(result.protein_id)

        self.output = BindingPocketPredictionOutput(proteins_with_binding_pockets=protein_ids)

    async def setup_jobs(self):
        self.jobs = []

        feature = SetupJobFeature()

        for protein_id in self.input.proteins_with_pdb:
            protein = Protein.objects.with_id(protein_id)

            job = await feature.handle(request=SetupJobRequest(
                experiment_id=self.experiment.id,
                protein_id=protein_id,
                job_name=f'Binding pocket prediction for protein {protein.name}'
            ))

            job = PocketPredictionJob.objects.with_id(job.job_id)

            self.jobs.append(job)

    async def prevalidate_jobs(self) -> List[JobValidationError]:
        validation_errors = []

        job: PocketPredictionJob
        for job in self.jobs:
            if not job.input_valid():
                validation_errors.append(
                    JobValidationError(
                        job_id=job.id,
                        msg=f'Job input is invalid. Setup job inputs manually'
                    )
                )

        return validation_errors

    @property
    def _input_parameter_type(self) -> Type[BindingPocketPredictionInput]:
        return BindingPocketPredictionInput

    @property
    def _output_parameter_type(self) -> Type[BindingPocketPredictionOutput]:
        return BindingPocketPredictionOutput
