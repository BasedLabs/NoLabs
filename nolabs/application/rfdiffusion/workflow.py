import uuid
from io import StringIO
from typing import List, Type, Optional, Dict, Any

from pydantic import BaseModel

from nolabs.application.rfdiffusion.worker_models import RunRfdiffusionRequest, RunRfdiffusionResponse
from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import Protein, JobId, JobName, ProteinName
from nolabs.domain.models.protein_design import RfdiffusionJob
from nolabs.infrastructure.mongo_connector import get_connection
from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler
from Bio.PDB import PDBParser, PDBIO


class RfDiffusionInput(BaseModel):
    contig: str
    timesteps: int
    proteins_with_pdb: List[uuid.UUID]
    number_of_designs: int = 1
    hotspots: str = ''
    inpaint: str = ''
    remove_chain: str = ''


class RfDiffusionOutput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class RfDiffusionComponent(Component[RfDiffusionInput, RfDiffusionOutput]):
    name = 'RFDiffusion Binder Design'
    description = 'Protein binder prediction using RFDiffusion (generates PDBs of binders)'

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
                          timesteps=inp.timesteps,
                          inpaint=inp.inpaint,
                          remove_chain=inp.remove_chain)
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
            inpaint=job.inpaint,
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
                if job.remove_chain:
                    pdb_content = self._remove_chain_from_pdb(pdb_content=pdb_content,
                                                              chain_id_to_remove=job.remove_chain)
                protein = Protein.create(
                    experiment=self.experiment_id,
                    name=ProteinName(f"{job.protein.name.value}-binder"),
                    pdb_content=pdb_content
                )
                protein.save()
                proteins.append(protein)
            job.set_result(binders=proteins)
            await job.save()

    def _remove_chain_from_pdb(self, pdb_content, chain_id_to_remove):
        """
        Remove a specified chain from PDB content provided as a string.

        Parameters:
        - pdb_content (str): The content of the PDB file as a string.
        - chain_id_to_remove (str): The ID of the chain to remove.

        Returns:
        - str: The modified PDB content as a string without the specified chain.
        """
        # Create a file-like object from the input string
        pdb_io = StringIO(pdb_content)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein_structure', pdb_io)

        for model in structure:
            # Collect chains to remove
            chains_to_remove = []
            for chain in model:
                # Check if the chain ID matches the one to remove
                if chain.id == chain_id_to_remove:
                    chains_to_remove.append(chain)
            for chain in chains_to_remove:
                model.detach_child(chain.id)

        output_io = StringIO()
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_io)
        modified_pdb_content = output_io.getvalue()

        return modified_pdb_content

    async def on_completion(
            self, inp: RfDiffusionInput, job_ids: List[uuid.UUID]
    ) -> Optional[RfDiffusionOutput]:
        protein_ids = []
        for job_id in job_ids:
            job: RfdiffusionJob = RfdiffusionJob.objects.with_id(job_id)
            for binder in job.binders:
                protein_ids.append(binder.id)
        return RfDiffusionOutput(proteins_with_pdb=protein_ids)


