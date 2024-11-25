__all__ = ["ProteinMPNNComponent"]

import uuid
from itertools import chain
from typing import List, Type, Optional, Dict, Any, Union

from pydantic import BaseModel

from nolabs.application.proteinmpnn.worker_models import (
    RunProteinMPNNPredictionInput,
    RunProteinMPNNPredictionOutput, )
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import JobId, JobName, Protein, ProteinName, ProteinId
from nolabs.domain.models.proteinmpnn import ProteinMPNNJob, ProteinMPNNResult
from nolabs.infrastructure.mongo_connector import get_connection
from nolabs.workflow.core.component import Component, TInput, TOutput
from nolabs.workflow.core.flow import ComponentFlowHandler


class ProteinMPNNComponentInput(BaseModel):
    num_seq_per_target: int = 2
    sampling_temp: float = 0.1
    seed: int = 37
    batch_size: int = 1
    proteins_with_pdb: Union[List[uuid.UUID], List[List[uuid.UUID]]]


class ProteinMPNNComponentOutput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]


class ProteinMPNNComponent(Component[ProteinMPNNComponentInput, ProteinMPNNComponentOutput]):
    name = "ProteinMPNN"
    description = "Protein design using ProteinMPNN. Predicts sequences from PDB backbones."

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
    async def on_finish(
            self, inp: ProteinMPNNComponentInput, job_ids: List[uuid.UUID]
    ) -> Optional[ProteinMPNNComponentOutput]:
        protein_ids = []

        for job_id in job_ids:
            job: ProteinMPNNJob = ProteinMPNNJob.objects.with_id(job_id)
            for result in job.results:
                protein_ids.append(result.protein.id)

        return ProteinMPNNComponentOutput(
            proteins_with_fasta=protein_ids
        )

    async def on_start(self, inp: ProteinMPNNComponentInput) -> List[uuid.UUID]:
        job_ids = []

        if not inp.proteins_with_pdb:
            return []

        if isinstance(inp.proteins_with_pdb[0], list):
            inp.proteins_with_pdb = list(chain.from_iterable(inp.proteins_with_pdb))

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
                job.set_input(protein=protein,
                              sampling_temp=inp.sampling_temp,
                              seed=inp.seed,
                              batch_size=inp.batch_size)
                await job.save()

            job_ids.append(job.id)

        return job_ids

    async def on_job_start(self, job_id: uuid.UUID):
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
            celery_task_name="design",
            celery_queue="proteinmpnn",
            input={"param": request.model_dump()},
        )

    async def on_job_finish(
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
            for item in job_result.fasta_contents:
                fasta_contents = ['>' + i.strip() for i in item.split('>') if i]
                for fasta_content in fasta_contents:
                    sequence = fasta_content.split('\n')[-1] # get sequence
                    upper_part = fasta_content[1:].split('\n')[0].split(', ') # get parts without sequence
                    d = {} # dictionary of upper parts
                    for comma_separated_part in upper_part:
                        if '=' in comma_separated_part:
                            name, value = comma_separated_part.split('=')
                            d[name] = value

                    protein_id = ProteinId(uuid.uuid4())
                    protein = job.protein.copy(id=protein_id)
                    protein.set_fasta(fasta_content=fasta_content)
                    protein.save(session=session)

                    sequence_result = ProteinMPNNResult(
                        id=uuid.uuid4(),
                        fasta_content=fasta_content,
                        sequence=sequence,
                        score= float(d.get("score")),
                        global_score=float(d.get("global_score")),
                        T=float(d.get("T")) if 'T' in d else None,
                        sample=float(d.get("sample")) if 'T' in d else None,
                        seq_recovery=float(d.get("seq_recovery")) if 'seq_recovery' in d else None,
                        protein=protein.id
                    )
                    sequence_result.save(session=session)
                    results.append(sequence_result)

            job.set_result(results=results)
            await job.save(session=session, cascade=True)
