import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.adaptyv_bio.protein_binding_screening_job import ProteinBindingScreeningJob
from nolabs.domain.models.common import JobId, JobName, Protein
from nolabs.infrastructure.settings import settings
from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler


class ProteinBindingScreeningInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]


class ProteinBindingScreeningOutput(BaseModel):
    pass


class ProteinBindingScreeningFlowHandler(ComponentFlowHandler):
    async def on_start(self, inp: ProteinBindingScreeningInput) -> List[uuid.UUID]:
        if not settings.adaptyv_bio_api_token:
            raise NoLabsException(ErrorCodes.adaptyv_bio_token_not_set)

        if not settings.adaptyv_bio_api_base:
            raise NoLabsException(ErrorCodes.adaptyv_bio_token_not_set)

        protein_ids = inp.proteins_with_fasta
        job: ProteinBindingScreeningJob = ProteinBindingScreeningJob.objects(
            proteins__in=protein_ids).first()
        if not job:
            # create a job
            job: ProteinBindingScreeningJob = ProteinBindingScreeningJob.create(
                id=JobId(uuid.uuid4()),
                name=JobName(f"{str(len(protein_ids))} proteins"),
                component=self.component_id
            )
            proteins = Protein.objects(id__in=protein_ids)
            job.set_proteins(proteins=proteins)
            await job.save()

        if len(job.proteins) == len(protein_ids):
            return [job.id]

        raise NoLabsException(ErrorCodes.proteins_are_part_of_another_job, "Proteins some or all of the proteins are used in another job/experiment")

    async def on_job_start(self, job_id: uuid.UUID):
        job: ProteinBindingScreeningJob = ProteinBindingScreeningJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        if not job.submitted:
            return self.cancel_job(job_id=job_id, reason="Open job and submit it manually")


class ProteinBindingScreeningComponent(
    Component[ProteinBindingScreeningInput, ProteinBindingScreeningOutput]):
    name = 'Protein binding screening'
    description = 'Protein binding screening (Adaptyv Bio)'

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return ProteinBindingScreeningInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return ProteinBindingScreeningOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlowHandler"]:
        return ProteinBindingScreeningFlowHandler
