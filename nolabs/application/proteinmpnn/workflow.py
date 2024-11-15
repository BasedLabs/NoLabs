__all__ = ["ProteinMPNNComponent"]

import uuid
from typing import List, Type, Optional, Dict, Any

from pydantic import BaseModel

from nolabs.application.proteinmpnn.worker_models import (
    RunProteinMPNNPredictionInput,
    RunProteinMPNNPredictionOutput,
)
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import JobId, JobName, Protein
from nolabs.domain.models.proteinmpnn import ProteinMPNNJob, ProteinMPNNResult
from nolabs.infrastructure.mongo_connector import get_connection
from nolabs.workflow.core.component import Component
from nolabs.workflow.core.flow import ComponentFlowHandler


class ProteinMPNNComponentInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]
    # Additional inputs if needed


class ProteinMPNNComponentOutput(BaseModel):
    generated_sequences: List[uuid.UUID]


class ProteinMPNNComponent(Component[ProteinMPNNComponentInput, ProteinMPNNComponentOutput]):
    name = "ProteinMPNN"
    description = "Protein design using ProteinMPNN"

    @property
    def input_parameter_type(self) -> Type[ProteinMPNNComponentInput]:
        return ProteinMPNNComponentInput

    @property
    def output_parameter_type(self) -> Type[ProteinMPNNComponentOutput]:
        return ProteinMPNNComponentOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlowHandler"]:
        return ProteinMPNNComponentFlowHandler


class ProteinMPNNComponentFlowHandler(ComponentFlowHandler):
    async def on_component_task(self, inp: ProteinMPNNComponentInput) -> List[uuid.UUID]:
        job_ids = []

        for protein_id in inp.proteins_with_pdb:
            protein: Protein = Protein.objects.with_id(protein_id)

            if not protein:
                raise NoLabsException(
                    ErrorCodes.protein_not_found, data={"protein_id": protein_id}
                )

            job = ProteinMPNNJob.objects(protein=protein.id).first()

            if not job:
                job = ProteinMPNNJob.create(
                    id=JobId(uuid.uuid4()),
                    name=JobName("ProteinMPNN design job"),
                    component=self.component_id,
                )
                job.set_input(protein=protein)
                await job.save()

            job_ids.append(job.id)

        return job_ids

    async def on_job_task(self, job_id: uuid.UUID):
        job: ProteinMPNNJob = ProteinMPNNJob.objects.with_id(job_id)

        input_errors = job.input_errors(throw=False)

        if input_errors:
            message = ", ".join(i.message for i in input_errors)
            raise NoLabsException(ErrorCodes.invalid_job_input, message=message)

        request = RunProteinMPNNPredictionInput(
            pdb_contents=job.protein.get_pdb(),
            is_homomer=job.is_homomer,
            chains_to_design=job.chains_to_design,
            fixed_positions=job.fixed_positions,
            num_seq_per_target=job.num_seq_per_target,
            sampling_temp=job.sampling_temp,
            seed=job.seed,
            batch_size=job.batch_size,
        )
        task_id = uuid.uuid4()
        job.set_task_id(task_id=str(task_id))
        await job.save()

        return await self.schedule(
            job_id=job.id,
            celery_task_name="run_protmpnn",
            celery_queue="proteinmpnn",
            input={"param": request.model_dump()},
        )

    async def on_job_completion(
        self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
    ):
        job: ProteinMPNNJob = ProteinMPNNJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        job_result = RunProteinMPNNPredictionOutput(**long_running_output)

        db = get_connection()
        session = db.client.start_session()
        with session.start_transaction():
            results = []
            for seq_data in job_result.sequences:
                sequence_result = ProteinMPNNResult(
                    id=uuid.uuid4(),
                    sequence=seq_data.sequence,
                    fasta_content=seq_data.fasta_content,
                    score=seq_data.score,
                    global_score=seq_data.global_score,
                    T=seq_data.T,
                    sample=seq_data.sample,
                    seq_recovery=seq_data.seq_recovery,
                )
                sequence_result.save(session=session)
                results.append(sequence_result)

            job.set_result(results=results)
            await job.save(session=session, cascade=True)
