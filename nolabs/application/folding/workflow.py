import uuid
from abc import ABC
from enum import Enum
from typing import Any, Dict, List, Optional, Type

from pydantic import BaseModel

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import JobId, JobName, Protein
from nolabs.domain.models.folding import FoldingJob
from nolabs.workflow.core.component import Component, TInput, TOutput
from nolabs.workflow.core.flow import ComponentFlowHandler


class FoldingComponentInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]


class FoldingComponentOutput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class FoldingBackendEnum(str, Enum):
    rosettafold = "rosettafold"
    esmfold = "esmfold"
    esmfold_light = "esmfold_light"


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


class FoldingComponentFlow(ComponentFlowHandler):
    job_timeout_seconds = 10.0
    component_timeout_seconds = 120.0

    backend: FoldingBackendEnum

    async def on_component_task(self, inp: FoldingComponentInput) -> List[uuid.UUID]:
        try:
            job_ids = []

            for protein_id in inp.proteins_with_fasta:
                protein = Protein.objects.with_id(protein_id)

                if not protein:
                    raise NoLabsException(ErrorCodes.protein_not_found)

                job: FoldingJob = FoldingJob.objects(protein=protein_id).first()

                if not job:
                    job_id = JobId(uuid.uuid4())
                    job_name = JobName(f"Folding of {protein.name}")
                    job = FoldingJob.create(
                        id=job_id,
                        name=job_name,
                        component=self.component_id,
                    )
                    job.set_inputs(protein=protein, backend=self.backend)
                    await job.save(cascade=True)

                job_ids.append(job.id)

            return [i for i in job_ids]
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.workflow_get_folding_jobs_failed) from e

    async def on_completion(
        self, inp: FoldingComponentInput, job_ids: List[uuid.UUID]
    ) -> FoldingComponentOutput:
        items = []

        for job_id in job_ids:
            job: FoldingJob = FoldingJob.objects.with_id(job_id)
            items.append(job.folded_protein.id)

        return FoldingComponentOutput(proteins_with_pdb=items)

    async def on_job_task(self, job_id: uuid.UUID):
        job: FoldingJob = FoldingJob.objects.with_id(job_id)

        input_errors = job.input_errors(throw=False)

        if input_errors:
            message = ", ".join(i.message for i in input_errors)

            raise NoLabsException(ErrorCodes.component_input_invalid, message=message)

        if not job.processing_required:
            return

        return await self.schedule(
            job_id=job.id,
            celery_task_name="inference",
            celery_queue="esmfold-light-service",
            input={"param": {"fasta_sequence": job.protein.get_fasta()}},
        )

    async def on_job_completion(
        self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
    ):
        job: FoldingJob = FoldingJob.objects.with_id(job_id)
        id = uuid.uuid4()
        folded_protein = job.protein.copy(id=id)
        folded_protein.set_pdb(long_running_output["pdb_content"])
        job.set_result(protein=folded_protein)
        await job.save(cascade=True)


class EsmfoldLightComponentFlow(FoldingComponentFlow):
    backend = FoldingBackendEnum.esmfold_light
