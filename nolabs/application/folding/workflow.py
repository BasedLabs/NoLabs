import uuid
from abc import ABC
from typing import List, Optional, Type

from domain.exceptions import ErrorCodes, NoLabsException
from prefect.states import Cancelled, Completed, Failed
from pydantic import BaseModel

from nolabs.application.folding.api_models import FoldingBackendEnum
from nolabs.domain.models.common import Experiment, JobId, JobName, Protein
from nolabs.domain.models.folding import FoldingJob
from nolabs.infrastructure.cel import cel as celery
from nolabs.microservices.esmfold_light.service.api_models import (
    InferenceInput)
from application.workflow import ComponentFlow
from application.workflow.component import Component, TInput, TOutput


class FoldingComponentInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]


class FoldingComponentOutput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class FoldingComponent(ABC, Component[FoldingComponentInput, FoldingComponentOutput]):
    name = "Folding"
    description = "Folding component"

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return FoldingComponentInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return FoldingComponentOutput


class EsmfoldComponent(FoldingComponent):
    name = "Esmfold"
    description = "Protein folding using Esmfold"

    @property
    def component_flow_type(self) -> Optional[Type["FoldingComponentFlow"]]:
        return FoldingComponentFlow


class EsmfoldLightComponent(FoldingComponent):
    name = "Esmfold light"
    description = "Protein folding using Esmfold light"

    @property
    def component_flow_type(self) -> Optional[Type["FoldingComponentFlow"]]:
        return EsmfoldLightComponentFlow


class RosettafoldComponent(FoldingComponent):
    name = "Rosettafold component"
    description = "Protein folding using Rosettafold"

    @property
    def backend(self) -> FoldingBackendEnum:
        return FoldingBackendEnum.rosettafold

    @property
    def component_flow_type(self) -> Optional[Type["FoldingComponentFlow"]]:
        return FoldingComponentFlow


class FoldingComponentFlow(ComponentFlow):
    job_timeout_seconds = 10.0
    component_timeout_seconds = 120.0

    backend: FoldingBackendEnum

    async def get_jobs(self, inp: FoldingComponentInput) -> List[uuid.UUID]:
        job_ids = []

        experiment = Experiment.objects.with_id(self.experiment_id)

        for protein_id in inp.proteins_with_fasta:
            protein = Protein.objects.with_id(protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            job: FoldingJob = FoldingJob.objects(protein=protein_id).first()

            if not job:
                job_id = JobId(uuid.uuid4())
                job_name = JobName(f"Folding of {protein.name}")
                job = FoldingJob(id=job_id, name=job_name, experiment=experiment)
                job.set_inputs(protein=protein, backend=self.backend)
                await job.save(cascade=True)

            job_ids.append(job.id)

        return [i for i in job_ids]

    async def gather_jobs(self, inp: FoldingComponentInput, job_ids: List[uuid.UUID]):
        items = []

        for job_id in job_ids:
            job: FoldingJob = FoldingJob.objects.with_id(job_id)
            items.append(job.protein.id)

        return FoldingComponentOutput(proteins_with_pdb=items)

    async def job_task(self, job_id: uuid.UUID):
        job: FoldingJob = FoldingJob.objects.with_id(job_id)

        input_errors = job.input_errors(throw=False)

        if input_errors:
            message = ", ".join(i.message for i in input_errors)

            return Cancelled(message=message)

        job_result = await celery.esmfold_light_inference(
            task_id=str(job.id),
            payload=InferenceInput(fasta_sequence=job.protein.get_fasta()),
        )

        protein = Protein.objects.with_id(job.protein.id)
        protein.set_pdb(job_result.pdb_content)
        protein.save()

        job.set_result(job.protein, job_result.pdb_content)

        await job.save()

        if not job.result_valid():
            return Failed(message="Result of job is invalid")

        return Completed()


class EsmfoldLightComponentFlow(FoldingComponentFlow):
    backend = FoldingBackendEnum.esmfold_light
