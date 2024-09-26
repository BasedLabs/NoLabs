__all__ = ["GetJobFeature", "RunJobFeature", "SetupJobFeature"]

from uuid import UUID, uuid4

from nolabs.application.blast.api_models import (
    HitModel,
    HspModel,
    JobResponse,
    JobResult,
    SetupJobRequest,
)
from nolabs.application.blast.services import BlastJobRunner
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.blast import BlastJob
from nolabs.domain.models.common import Experiment, JobId, JobName, Protein


def map_job_to_response(job: BlastJob) -> JobResponse:
    return JobResponse(
        job_id=job.id,
        job_name=job.name.value,
        protein_id=job.protein.id,
        descriptions=job.descriptions,
        alignments=job.alignments,
        hitlist_size=job.hitlist_size,
        expect=job.expect,
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
                                midline=hsp.midline,
                            )
                            for hsp in hit.hsps
                        ],
                    )
                    for hit in res.hits
                ],
            )
            for res in job.result
        ],
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

        job_id = JobId(request.job_id if request.job_id else uuid4())
        job_name = JobName(request.job_name if request.job_name else "New BLAST job")

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job: BlastJob = BlastJob.objects.with_id(job_id.value)

        if not job:
            job = BlastJob(id=job_id, name=job_name, experiment=experiment)

        protein = Protein.objects.with_id(request.protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        job.set_input(
            protein=protein,
            job_type="blastp",
            descriptions=request.descriptions,
            alignments=request.alignments,
            hitlist_size=request.hitlist_size,
            expect=request.expect,
        )

        await job.save(cascade=True)

        return map_job_to_response(job)


class RunJobFeature:
    """
    Use case - start job.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job: BlastJob = BlastJob.objects.get(id=job_id)

        job_runner = BlastJobRunner()
        await job_runner.run(job=job)

        return map_job_to_response(job)
