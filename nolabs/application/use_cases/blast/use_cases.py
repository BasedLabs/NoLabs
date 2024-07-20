__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature',
    'GetJobStatusFeature'
]

from typing import List
from uuid import UUID

import blast_query_microservice
from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.blast.api_models import JobResponse, JobResult, \
    SetupJobRequest, GetJobStatusResponse
from nolabs.domain.models.common import JobId, Experiment, JobName, Protein, Ligand
from nolabs.domain.models.blast import BlastJob, BlastJobResult
from nolabs.utils import generate_uuid


def map_job_to_response(job: BlastJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        samples_per_complex=job.samples_per_complex,
        protein_id=job.protein.iid.value,
        result=[
            JobResult(
                complex_id=res.complex_id,
                sdf_content=res.sdf_content.decode('utf-8'),
                minimized_affinity=res.minimized_affinity,
                scored_affinity=res.scored_affinity,
                confidence=res.confidence
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

        job.set_input(protein=protein)
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

        result: List[BlastJobResult] = []

        request = blast_query_microservice.blast(pdb_contents=job.protein.get_pdb(),
                                                                     samples_per_complex=job.samples_per_complex,
                                                                     job_id=str(job_id.value),
                                                                     )

        try:
            job.started()
            await job.save()
            response = self._blast.predict_run_docking_post(
                run_diff_dock_prediction_request=request, _request_timeout=(1000.0, 1000.0))

            if not response.success:
                raise NoLabsException(ErrorCodes.blast_api_error, response.message)

        finally:
            job.finished()
            await job.save()

        return map_job_to_response(job)
