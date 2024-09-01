import uuid
from abc import ABC, abstractmethod
from typing import List, Type, Optional, Dict, Any

from celery.result import AsyncResult
from pydantic import BaseModel

from nolabs.application.use_cases.folding.api_models import FoldingBackendEnum
from nolabs.application.workflow import ComponentTask
from nolabs.application.workflow.component import Component, JobValidationError, TOutput, TInput
from nolabs.domain.models.common import Protein, JobId, JobName, Experiment
from nolabs.domain.models.folding import FoldingJob
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.infrastructure.celery_tasks import esmfold_light_inference
from nolabs.infrastructure.celery_worker import celery_app, send_selery_task
from nolabs.microservices.esmfold_light.service.api_models import InferenceInput, InferenceOutput


class FoldingComponentInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]
    backend: Optional[FoldingBackendEnum] = FoldingBackendEnum.esmfold_light


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
    def component_task_type(self) -> Optional[Type['FoldingComponentTask']]:
        return FoldingComponentTask


class EsmfoldLightComponent(FoldingComponent):
    name = 'Esmfold light'
    description = 'Protein folding using Esmfold light'

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.esmfold_light

    @property
    def component_task_type(self) -> Optional[Type['FoldingComponentTask']]:
        return FoldingComponentTask


class RosettafoldComponent(FoldingComponent):
    name = 'Rosettafold component'
    description = 'Protein folding using Rosettafold'

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.rosettafold

    @property
    def component_task_type(self) -> Optional[Type['FoldingComponentTask']]:
        return FoldingComponentTask


class FoldingComponentTask(ComponentTask):
    job_timeout_seconds = 10.0
    component_timeout_seconds = 120.0

    async def get_jobs(self, inp: FoldingComponentInput) -> List[uuid.UUID]:
        job_ids = []

        experiment = Experiment.objects.with_id(self.extra['experiment_id'])

        for protein_id in inp.proteins_with_fasta:
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

            job.set_inputs(protein=protein, backend=inp.backend)
            await job.save(cascade=True)

            job_ids.append(job.id)

        return [i for i in job_ids]

    async def post_execute(self, inp: FoldingComponentInput, job_ids: List[uuid.UUID]):
        items = []

        for job_id in job_ids:
            job: FoldingJob = FoldingJob.objects.with_id(job_id)
            items.append(job.protein.id)

        return FoldingComponentOutput(
            proteins_with_pdb=items
        )

    async def execute_job_task(self, job_id: uuid.UUID):
        job: FoldingJob = FoldingJob.objects.with_id(job_id)

        job.input_errors(throw=True)

        async_result: AsyncResult = send_selery_task(name=esmfold_light_inference,
                                                     payload=InferenceInput(fasta_sequence=job.protein.get_fasta()))
        job_result: Dict[str, Any] = await self.celery_wait_async(async_result)
        job_result: InferenceOutput = InferenceOutput(**job_result)

        job.set_result(job.protein, job_result.pdb_content)

        await job.save()
