__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature',
    'GetJobStatusFeature'
]

from uuid import UUID

import p2rank_microservice
from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.binding_pockets.api_models import GetJobStatusResponse, JobResponse, \
    SetupJobRequest
from nolabs.refined.domain.models.common import Protein, JobId, JobName, Experiment
from nolabs.refined.domain.models.pocket_prediction import PocketPredictionJob
from nolabs.utils import generate_uuid


def map_job_to_response(job: PocketPredictionJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        protein_id=job.protein.iid.value,
        result=job.pocket_ids
    )


class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        job_id = JobId(job_id)
        job: PocketPredictionJob = PocketPredictionJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return map_job_to_response(job)


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """

    async def handle(self, request: SetupJobRequest) -> JobResponse:
        assert request

        job_id = JobId(request.job_id or generate_uuid())
        job_name = JobName(request.job_name or 'New binding pocket prediction job')

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job: PocketPredictionJob = PocketPredictionJob.objects(Q(id=job_id.value) | Q(name=job_name.value)).first()

        if not job:
            job = PocketPredictionJob(
                id=job_id,
                name=job_name,
                experiment=experiment
            )

        protein = Protein.objects.with_id(request.protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        job.set_input(protein=protein)
        job.save(cascade=True)

        return map_job_to_response(job)


class RunJobFeature:
    """
    Use case - start job.
    """
    _api: p2rank_microservice.DefaultApi

    def __init__(self, api: p2rank_microservice.DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: PocketPredictionJob = PocketPredictionJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        if not job.protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty,
                                  'Cannot run binding pockets prediction on an empty protein pdb')

        protein = job.protein

        response = self._api.predict_run_p2rank_post(
            run_p2_rank_prediction_request=p2rank_microservice.RunP2RankPredictionRequest(
                job_id=str(job_id),
                pdb_contents=protein.get_pdb()
            )
        )

        job.set_result(protein=protein, pocket_ids=response.pocket_ids)
        job.save(cascade=True)

        protein.set_binding_pockets(response.pocket_ids)
        protein.save()

        return map_job_to_response(job)


class GetJobStatusFeature:
    _api: p2rank_microservice.DefaultApi

    def __init__(self, api: p2rank_microservice.DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        job: PocketPredictionJob = PocketPredictionJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        result = self._api.is_job_running_job_job_id_is_running_get(job_id=str(job_id))

        return GetJobStatusResponse(
            running=result.is_running,
            result_valid=job.result_valid()
        )
