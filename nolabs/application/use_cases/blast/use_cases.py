__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature',
    'GetJobStatusFeature'
]

from typing import List
from uuid import UUID

import blast_query_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.blast.api_models import JobResponse, JobResult, \
    SetupJobRequest, GetJobStatusResponse, HitModel, HspModel
from nolabs.domain.models.common import JobId, Experiment, JobName, Protein
from nolabs.domain.models.blast import BlastJob, Hit, Hsp, BlastJobResult
from nolabs.utils import generate_uuid


def map_job_to_response(job: BlastJob) -> JobResponse:
    return JobResponse(
        job_id=job.id,
        job_name=job.name.value,
        protein_id=job.protein.id,
        result=[
            JobResult(
                protein_id=res.protein_id,
                program=res.program,
                database=res.database,
                query_id=res.query_id,
                query_def=res.query_def,
                query_len=res.query_len,
                hits=[
                    HitModel(
                        num=hit.num,
                        id=hit.id,
                        definition=hit.definition,
                        accession=hit.accession,
                        length=hit.length,
                        hsps=[
                            HspModel(
                                num=hsp.num,
                                bit_score=hsp.bit_score,
                                score=hsp.score,
                                evalue=hsp.evalue,
                                query_from=hsp.query_from,
                                query_to=hsp.query_to,
                                hit_from=hsp.hit_from,
                                hit_to=hsp.hit_to,
                                query_frame=hsp.query_frame,
                                hit_frame=hsp.hit_frame,
                                identity=hsp.identity,
                                positive=hsp.positive,
                                gaps=hsp.gaps,
                                align_len=hsp.align_len,
                                qseq=hsp.qseq,
                                hseq=hsp.hseq,
                                midline=hsp.midline
                            )
                            for hsp in hit.hsps
                        ]
                    )
                    for hit in res.hits
                ]
            )
            for res in job.result
        ]
    )



class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        job_id = JobId(job_id)
        job: BlastJob = BlastJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return map_job_to_response(job)


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """

    async def handle(self, request: SetupJobRequest) -> JobResponse:
        assert request

        job_id = JobId(request.job_id if request.job_id else generate_uuid())
        job_name = JobName(request.job_name if request.job_name else 'New protein ligand DIFFDOCK binding job')

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job: BlastJob = BlastJob.objects.with_id(job_id.value)

        if not job:
            job = BlastJob(
                id=job_id,
                name=job_name,
                experiment=experiment
            )

        protein = Protein.objects.with_id(request.protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        job.set_input(protein=protein,
                      job_type="blastp",
                      descriptions=request.descriptions,
                      alignments=request.alignments,
                      hitlist_size=request.hitlist_size,
                      expect=request.expect)

        await job.save(cascade=True)

        return map_job_to_response(job)


class GetJobStatusFeature:
    """
    Use case - set job status.
    """
    _blast = blast_query_microservice.DefaultApi

    def __init__(self,
                 blast: blast_query_microservice.DefaultApi):
        self._blast = blast

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        job: BlastJob = BlastJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        response = self._blast.is_job_running_job_job_id_is_running_get(
            job_id=str(job.iid.value)
        )
        return GetJobStatusResponse(
            running=response.is_running,
            result_valid=job.result_valid()
        )


class RunJobFeature:
    """
    Use case - start job.
    """
    _blast = blast_query_microservice.DefaultApi

    def __init__(self, blast: blast_query_microservice.DefaultApi):
        self._blast = blast

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: BlastJob = BlastJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        request = blast_query_microservice.SequenceQuery(
            sequence=job.protein.get_amino_acid_sequence(),
            type=blast_query_microservice.BlastType.BLASTP,
            descriptions=job.descriptions,
            alignments=job.alignments,
            hitlist_size=job.hitlist_size,
            expect=job.expect,
            job_id=str(job.id)
        )

        try:
            job.started()
            await job.save()
            blast_output = self._blast.blast_blast_post(
                sequence_query=request, _request_timeout=(1000.0, 1000.0))

            if not blast_output:
                raise NoLabsException(ErrorCodes.blast_api_error)

            # Process response to create BlastJobResult
            hits = create_hits_and_hsps(blast_output)

            blast_job_result = BlastJobResult.create(
                protein_id=job.protein.id,
                program=blast_output['BlastOutput']['BlastOutput_program'],
                database=blast_output['BlastOutput']['BlastOutput_db'],
                query_id=blast_output['BlastOutput']['BlastOutput_query-ID'],
                query_def=blast_output['BlastOutput']['BlastOutput_query-def'],
                query_len=int(blast_output['BlastOutput']['BlastOutput_query-len']),
                hits=hits
            )

            job.set_result([blast_job_result])

        finally:
            job.finished()
            await job.save()

        return map_job_to_response(job)

        # Function to create hits and hsps from BLAST output
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
