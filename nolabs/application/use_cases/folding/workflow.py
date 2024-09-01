import uuid
from abc import ABC, abstractmethod
from typing import List, Type, Optional, Dict, Any

from celery.result import AsyncResult
from prefect import task
from pydantic import BaseModel

from nolabs.application.use_cases.folding.api_models import FoldingBackendEnum
from nolabs.application.workflow import SetupTask, OutputTask, ExecuteJobTask
from nolabs.application.workflow.component import Component, JobValidationError, TOutput, TInput
from nolabs.domain.models.common import Protein, JobId, JobName, Experiment
from nolabs.domain.models.folding import FoldingJob
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.infrastructure.celery_tasks import esmfold_light_inference
from nolabs.infrastructure.celery_worker import celery_app, send_selery_task
from nolabs.microservices.esmfold_light.service.api_models import InferenceInput, InferenceOutput


class FoldingComponentInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]


class FoldingComponentOutput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class FoldingComponent(ABC, Component[FoldingComponentInput, FoldingComponentOutput]):
    name = 'Folding'
    description = 'Folding component'

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return FoldingComponentInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return FoldingComponentOutput

    @property
    def output_task_type(self) -> Type['OutputTask']:
        return GatherJobOutputsTask

    @property
    def job_task_type(self) -> Optional[Type['ExecuteJobTask']]:
        return ExecuteFoldingJobTask

    @property
    def setup_task_type(self) -> Type['SetupTask']:
        return SetupFoldingJobsTask

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


class EsmfoldComponent(FoldingComponent):
    name = 'Esmfold'
    description = 'Protein folding using Esmfold'

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.esmfold

    @property
    def job_task_type(self) -> Optional[Type['ExecuteJobTask']]:
        return ExecuteFoldingJobTask


class EsmfoldLightComponent(FoldingComponent):
    name = 'Esmfold light'
    description = 'Protein folding using Esmfold light'

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.esmfold_light

    @property
    def job_operator_type(self) -> Optional[Type['ExecuteJobTask']]:
        return ExecuteFoldingJobTask


class RosettafoldComponent(FoldingComponent):
    name = 'Rosettafold component'
    description = 'Protein folding using Rosettafold'

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.rosettafold

    @property
    def job_operator_type(self) -> Optional[Type['ExecuteJobTask']]:
        return ExecuteFoldingJobTask


class SetupFoldingJobsTask(SetupTask):
    timeout_seconds = 10.0

    async def execute(self) -> List[uuid.UUID]:
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

            job_ids.append(job.id)

        component.job_ids = job_ids
        component.save()

        return [i for i in job_ids]


class ExecuteFoldingJobTask(ExecuteJobTask):
    timeout_seconds = 10.0

    async def execute(self, job_id: uuid.UUID) -> Optional[BaseModel]:
        job: FoldingJob = FoldingJob.objects.with_id(job_id)

        job.input_errors(throw=True)

        async_result: AsyncResult = send_selery_task(name=esmfold_light_inference,
                                                     payload=InferenceInput(fasta_sequence=job.protein.get_fasta()))
        job_result: Dict[str, Any] = await self.celery_wait_async(async_result)
        job_result: InferenceOutput = InferenceOutput(**job_result)

        job.set_result(job.protein, job_result.pdb_content)

        await job.save()

        return


class GatherJobOutputsTask(OutputTask):
    async def execute(self) -> Optional[BaseModel]:
        component = self.get_component()

        items = []

        for job_id in component.job_ids:
            job: FoldingJob = FoldingJob.objects.with_id(job_id)
            items.append(job.protein.id)

        self.setup_output(FoldingComponentOutput(
            proteins_with_pdb=items
        ))
        return None
