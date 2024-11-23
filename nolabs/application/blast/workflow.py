import uuid
from typing import List, Type, Optional, Dict, Any

from pydantic import BaseModel

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.blast import BlastJob, Hsp, Hit, BlastJobResult
from nolabs.domain.models.common import Protein, JobId, JobName
from nolabs.workflow.core.component import Component, TOutput, TInput
from nolabs.workflow.core.flow import ComponentFlowHandler


class BlastComponentInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]
    descriptions: int = 10
    alignments: int = 10
    hitlist_size: int = 10
    expect: int = 10


class BlastComponentOutput(BaseModel):
    ...


def create_hits_and_hsps(blast_output: dict) -> List[Hit]:
    hits = []
    iteration_hits = blast_output['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
    for hit in iteration_hits:
        hsps = []
        hsp=hit['Hit_hsps']['Hsp']
        hsps.append(Hsp(
            num=int(hsp['Hsp_num']),
            bit_score=float(hsp['Hsp_bit-score']),
            score=int(hsp['Hsp_score']),
            evalue=float(hsp['Hsp_evalue']),
            query_from=int(hsp['Hsp_query-from']),
            query_to=int(hsp['Hsp_query-to']),
            hit_from=int(hsp['Hsp_hit-from']),
            hit_to=int(hsp['Hsp_hit-to']),
            query_frame=int(hsp['Hsp_query-frame']),
            hit_frame=int(hsp['Hsp_hit-frame']),
            identity=int(hsp['Hsp_identity']),
            positive=int(hsp['Hsp_positive']),
            gaps=int(hsp['Hsp_gaps']),
            align_len=int(hsp['Hsp_align-len']),
            qseq=hsp['Hsp_qseq'],
            hseq=hsp['Hsp_hseq'],
            midline=hsp['Hsp_midline']
        ))
        hits.append(Hit(
            num=int(hit['Hit_num']),
            id=hit['Hit_id'],
            definition=hit['Hit_def'],
            accession=hit['Hit_accession'],
            length=int(hit['Hit_len']),
            hsps=hsps
        ))
    return hits


class BlastComponentFlowHandler(ComponentFlowHandler):
    async def on_component_task(self, inp: BlastComponentInput) -> List[uuid.UUID]:
        jobs = []
        for protein_id in inp.proteins_with_fasta:
            job: BlastJob = BlastJob.objects(protein=protein_id).first()
            if not job:
                protein = Protein.objects.with_id(protein_id)
                job = BlastJob.create(
                    id=JobId(uuid.uuid4()),
                    name=JobName(f'Blast job for protein {protein.name}'),
                    component=self.component_id
                )
                job.set_input(
                    protein=protein,
                    job_type="blastp",
                    descriptions=inp.descriptions,
                    alignments=inp.alignments,
                    hitlist_size=inp.hitlist_size,
                    expect=inp.expect
                )
                await job.save()
            jobs.append(job.id)
        return jobs

    async def on_job_task(self, job_id: uuid.UUID):
        job: BlastJob = BlastJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return await self.schedule(
            job_id=job_id,
            celery_task_name="run_blast",
            celery_queue="blast",
            input={
                'param': {
                    'query': job.protein.get_amino_acid_sequence(),
                    'descriptions': job.descriptions,
                    'alignments': job.alignments,
                    'hitlist_size': job.hitlist_size,
                    'expect': job.expect
                }
            }
        )

    async def on_job_completion(
            self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
    ):
        job: BlastJob = BlastJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        hits = create_hits_and_hsps(long_running_output)

        blast_job_result = BlastJobResult.create(
            protein_id=job.protein.id,
            program=long_running_output['BlastOutput']['BlastOutput_program'],
            database=long_running_output['BlastOutput']['BlastOutput_db'],
            query_id=long_running_output['BlastOutput']['BlastOutput_query-ID'],
            query_def=long_running_output['BlastOutput']['BlastOutput_query-def'],
            query_len=int(long_running_output['BlastOutput']['BlastOutput_query-len']),
            hits=hits
        )

        job.set_result([blast_job_result])
        await job.save()


class BlastComponent(Component[BlastComponentInput, BlastComponentOutput]):
    name = 'Blast'
    description = 'Finds similar sequences for proteins and nucleotides.'

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return BlastComponentInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return BlastComponentOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlowHandler"]:
        return BlastComponentFlowHandler
