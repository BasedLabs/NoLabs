import uuid
from typing import List, Type, Optional, Dict, Any

from pydantic import BaseModel

from nolabs.application.rfdiffusion.worker_models import RunRfdiffusionRequest, RunRfdiffusionResponse
from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import Protein, JobId, JobName, ProteinName
from nolabs.domain.models.protein_design import RfdiffusionJob
from nolabs.infrastructure.mongo_connector import get_connection
from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler


class RfDiffusionInput(BaseModel):
    contig: str
    timesteps: int
    proteins_with_pdb: List[uuid.UUID]
    number_of_designs: int = 1
    hotspots: Optional[str] = None


class RfDiffusionOutput(BaseModel):
    binder_proteins: List[uuid.UUID]


class RfDiffusionComponent(Component[RfDiffusionInput, RfDiffusionOutput]):
    name = 'Protein binder design'
    description = 'Protein binder prediction using Rfdiffusion'

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return RfDiffusionInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return RfDiffusionOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlowHandler"]:
        return RfDiffusionFlowHandler

class RfDiffusionFlowHandler(ComponentFlowHandler):
    async def on_component_task(self, inp: RfDiffusionInput) -> List[uuid.UUID]:
        job_ids = []

        for protein_id in inp.proteins_with_pdb:
            protein = Protein.objects.with_id(protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            job: RfdiffusionJob = RfdiffusionJob.objects(protein=protein_id).first()

            if not job:
                job_id = JobId(uuid.uuid4())
                job_name = JobName(f"Folding of {protein.name}")
                job = RfdiffusionJob.create(
                    id=job_id,
                    name=job_name,
                    component=self.component_id,
                )
                job.set_input(protein=protein,
                              contig=inp.contig,
                              number_of_designs=inp.number_of_designs,
                              hotspots=inp.hotspots,
                              timesteps=inp.timesteps)
            job.set_input(protein=protein,
                          contig=inp.contig,
                          number_of_designs=inp.number_of_designs,
                          hotspots=inp.hotspots,
                          timesteps=inp.timesteps)
            await job.save(cascade=True)

            job_ids.append(job.id)

        return [i for i in job_ids]

    async def on_job_task(self, job_id: uuid.UUID):
        job: RfdiffusionJob = RfdiffusionJob.objects.with_id(job_id)

        input_errors = job.input_errors(throw=False)

        if input_errors:
            message = ", ".join(i.message for i in input_errors)

            raise NoLabsException(ErrorCodes.component_input_invalid, message=message)

        input = RunRfdiffusionRequest(
            pdb_content=job.protein.get_pdb(),
            contig=job.contig,
            hotspots=job.hotspots,
            timesteps=job.timesteps,
            number_of_designs=job.number_of_designs
        )

        return await self.schedule(job_id=job_id,
                                   celery_task_name="design",
                                   celery_queue="rfdiffusion",
                                   input={'param': input.model_dump()}) # TODO move names of tasks and queues, remove param (single input)


    async def on_job_completion(
            self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
    ):
        job: RfdiffusionJob = RfdiffusionJob.objects.with_id(job_id)
        output = RunRfdiffusionResponse(**long_running_output)
        if output.errors and not output.pdbs_content:
            raise NoLabsException(ErrorCodes.job_execution_failed, ", ".join(output.errors))
        proteins = []
        db = get_connection()
        session = db.client.start_session()
        with session.start_transaction():
            for pdb_content in output.pdbs_content:
                protein = Protein.create(
                    experiment=self.experiment_id,
                    name=ProteinName(f"{job.protein.name.value}-binder"),
                    pdb_content=pdb_content
                )
                protein.save()
                proteins.append(protein)
            job.set_result(binders=proteins)
            await job.save()

    async def on_completion(
            self, inp: RfDiffusionInput, job_ids: List[uuid.UUID]
    ) -> Optional[RfDiffusionOutput]:
        protein_ids = []
        for job_id in job_ids:
            job: RfdiffusionJob = RfdiffusionJob.objects.with_id(job_id)
            for binder in job.binders:
                protein_ids.append(binder.id)
        return RfDiffusionOutput(binder_proteins=protein_ids)


