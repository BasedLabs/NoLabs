import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.folding.api_models import FoldingBackendEnum, SetupJobRequest
from nolabs.refined.application.use_cases.folding.use_cases import SetupJobFeature, RunJobFeature, GetJobFeature
from nolabs.refined.domain.models.common import ExperimentId, Experiment, Protein
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import PythonComponent


class FoldingComponentInput(BaseModel):
    protein_ids: List[uuid.UUID]
    experiment_id: uuid.UUID


class FoldingComponentOutput(BaseModel):
    protein_ids: List[uuid.UUID]


class FoldingComponent(PythonComponent[FoldingComponentInput, FoldingComponentOutput]):
    name = 'Folding'

    async def _execute(self):
        setup_job_feature = SetupJobFeature()
        run_job_feature = RunJobFeature(
            esmfold=InfrastructureDependencies.esmfold_microservice()
        )
        get_job_feature = GetJobFeature()

        experiment_id = ExperimentId(self.input.experiment_id)

        experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        for protein_id in self.input.protein_ids:
            protein = Protein.objects.with_id(protein_id)

            protein.save()

            result = await setup_job_feature.handle(request=SetupJobRequest(
                experiment_id=experiment_id.value,
                backend=FoldingBackendEnum.esmfold,
                proteins=[protein.id],
                job_id=None,
                job_name=f'Folding {protein.name.fasta_name}'
            ))

            self._append_job_id(result.job_id)

        for job_id in self.job_ids:
            await run_job_feature.handle(job_id=job_id)

        items: List[uuid.UUID] = []

        for job_id in self.job_ids:
            get_result = await get_job_feature.handle(job_id=job_id)

            for protein_id in get_result.proteins:
                items.append(protein_id)

        self.output = FoldingComponentOutput(
            protein_ids=items
        )

    async def restore_output(self):
        get_job_feature = GetJobFeature()

        job_ids = self.get_job_ids()

        self.job_ids = job_ids

        items: List[uuid.UUID] = []

        for job_id in self.job_ids:
            get_result = await get_job_feature.handle(job_id=job_id)

            for protein_id in get_result.proteins:
                items.append(protein_id)

        self.output = FoldingComponentOutput(
            protein_ids=items
        )

    @property
    def _input_parameter_type(self) -> Type[FoldingComponentInput]:
        return FoldingComponentInput

    @property
    def _output_parameter_type(self) -> Type[FoldingComponentOutput]:
        return FoldingComponentOutput
