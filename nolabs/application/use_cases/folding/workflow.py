import uuid
from abc import ABC, abstractmethod
from typing import List, Type, Optional, Any

from airflow.models import BaseOperator
from airflow.utils.context import Context
from pydantic import BaseModel

from nolabs.application.use_cases.folding.api_models import FoldingBackendEnum
from nolabs.application.workflow import SetupOperator, OutputOperator
from nolabs.application.workflow.component import Component, JobValidationError, TOutput, TInput
from nolabs.domain.models.common import Protein, JobId, JobName, Experiment
from nolabs.domain.models.folding import FoldingJob
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.job_services.esmfold_light.operators.api_models import PredictFoldingJobRequest, PredictFoldingJobResponse
from nolabs.job_services.esmfold_light.operators.operator import RunEsmFoldLightOperator


class FoldingComponentInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]


class FoldingComponentOutput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class FoldingComponent(ABC, Component[FoldingComponentInput, FoldingComponentOutput]):
    @property
    def input_parameter_type(self) -> Type[TInput]:
        return FoldingComponentInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return FoldingComponentOutput

    @property
    def setup_operator_type(self) -> Type[BaseOperator]:
        return SetupFoldingJobsOperator

    @property
    def output_operator_type(self) -> Type[BaseOperator]:
        return GatherJobOutputsOperator

    @property
    @abstractmethod
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.esmfold_light

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


class EsmfoldComponent(FoldingComponent[FoldingComponentInput, FoldingComponentOutput]):
    name = 'Esmfold'
    description = 'Protein folding using Esmfold'

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.esmfold

    @property
    def job_operator_type(self) -> Optional[Type[BaseOperator]]:
        return RunEsmFoldLightOperator


class EsmfoldLightComponent(FoldingComponent):
    name = 'Esmfold light'
    description = 'Protein folding using Esmfold light'

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.esmfold_light

    @property
    def job_operator_type(self) -> Optional[Type[BaseOperator]]:
        return RunEsmFoldLightOperator


class RosettafoldComponent(FoldingComponent):
    name = 'Rosettafold component'
    description = 'Protein folding using Rosettafold'

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.rosettafold


class SetupFoldingJobsOperator(SetupOperator):
    async def execute_async(self, context: Context) -> List[str]:
        job_ids = []

        component: FoldingComponent = FoldingComponent.get(id=self.component_id)

        experiment = Experiment.objects.with_id(self.extra['experiment_id'])

        for protein_id in component.input_value.proteins_with_fasta:
            protein = Protein.objects.with_id(protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            job_id = JobId(uuid.uuid4())
            job_name = JobName('New folding job')
            job = FoldingJob(
                id=job_id,
                name=job_name,
                experiment=experiment
            )

            job.set_inputs(protein=protein, backend=component.backend)
            await job.save(cascade=True)

            self.set_job_input(job.id, PredictFoldingJobRequest(protein_sequence=protein.get_amino_acid_sequence()))

        component.job_ids = job_ids
        component.save()

        return [str(i) for i in job_ids]


class GatherJobOutputsOperator(OutputOperator):
    async def execute_async(self, context: Context) -> Any:
        component = self.get_component()

        items = []

        for job_id in component.job_ids:
            job: FoldingJob = FoldingJob.objects.with_id(job_id)
            job_result = PredictFoldingJobResponse(**self.get_job_output(job.id))
            job.set_result(job.protein, job_result.pdb_content)

            items.append(job.protein.iid.value)

        self.setup_output(FoldingComponentOutput(
            proteins_with_pdb=items
        ))
