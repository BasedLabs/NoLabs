import uuid
from typing import List, Type

from bson import ObjectId
from pydantic import BaseModel

from nolabs.application.adaptyv_bio.api_proxy import AdaptyvBioProteinAffinityCharacterizationApi
from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.adaptyv_bio.protein_affinity_characterization_job import ProteinAffinityCharacterizationJob
from nolabs.domain.models.common import JobId, JobName, Protein
from nolabs.infrastructure.settings import settings
from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler


class ProteinAffinityCharacterizationInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]


class ProteinAffinityCharacterizationOutput(BaseModel):
    pass


class ProteinAffinityCharacterizationFlowHandler(ComponentFlowHandler):
    async def on_start(self, inp: ProteinAffinityCharacterizationInput) -> List[uuid.UUID]:
        if not settings.adaptyv_bio_api_token:
            raise NoLabsException(ErrorCodes.adaptyv_bio_token_not_set)

        if not settings.adaptyv_bio_api_base:
            raise NoLabsException(ErrorCodes.adaptyv_bio_token_not_set)

        protein_ids = inp.proteins_with_fasta
        job: ProteinAffinityCharacterizationJob = ProteinAffinityCharacterizationJob.objects(
            proteins__in=protein_ids).first()
        if not job:
            # create a job
            job: ProteinAffinityCharacterizationJob = ProteinAffinityCharacterizationJob.create(
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
        job: ProteinAffinityCharacterizationJob = ProteinAffinityCharacterizationJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        if not job.submitted:
            raise NoLabsException(ErrorCodes.adaptyv_bio_job_was_not_submitted, "Open job and submit it manually")


class ProteinAffinityCharacterizationComponent(
    Component[ProteinAffinityCharacterizationInput, ProteinAffinityCharacterizationOutput]):
    name = 'Protein affinity characterization'
    description = 'Protein affinity characterization (Adaptyv bio)'

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return ProteinAffinityCharacterizationInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return ProteinAffinityCharacterizationOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlowHandler"]:
        return ProteinAffinityCharacterizationFlowHandler
