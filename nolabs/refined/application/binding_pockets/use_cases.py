__all__ = [
    'SetManualBindingPocketsFeature',
    'PredictBindingPocketsFeature',
    'GetBindingPocketsFeature'
]

from typing import List
from uuid import UUID

import p2rank_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.binding_pockets.api_models import BindingPocketsResponse, \
    PredictBindingPocketsRequest, GetJobStatusResponse, JobResponse, SetupJobRequest
from nolabs.refined.domain.models import Protein, JobId, JobName, Experiment
from nolabs.refined.domain.models.pocket_prediction import PocketPredictionJob
from nolabs.utils import generate_uuid


def map_job_to_response(job: PocketPredictionJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        protein_id=job.protein.iid.value,
        binding_pockets=job.pocket_ids
    )


class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            job_id = JobId(job_id)
            job: PocketPredictionJob = PocketPredictionJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            return map_job_to_response(job)
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e

            raise NoLabsException(ErrorCodes.unknown_exception) from e


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """
    async def handle(self, request: SetupJobRequest) -> JobResponse:
        try:
            assert request

            job_id = JobId(request.job_id if request.job_id else generate_uuid())
            job_name = JobName(request.job_name if request.job_name else 'New binding pocket prediction job')

            experiment = Experiment.objects.with_id(id=request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job: PocketPredictionJob = PocketPredictionJob.objects.with_id(job_id.value)

            if not job:
                job = PocketPredictionJob(
                    id=job_id,
                    name=job_name,
                    experiment=experiment
                )

            protein = Protein.objects.get(id=request.protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            job.set_input(protein=protein)
            job.save(cascade=True)

            return map_job_to_response(job)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e




class SetManualBindingPocketsFeature:
    async def handle(self, protein_id: UUID, binding_pockets: List[int]):
        protein = Protein.objects.with_id(id=protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        protein.set_manual_binding_pockets(binding_pockets=binding_pockets)
        protein.save()


class GetJobStatusFeature:
    _api: p2rank_microservice.DefaultApi

    def __init__(self, api: p2rank_microservice.DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        job: PocketPredictionJob = PocketPredictionJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        result = self._api.is_job_running_job_job_id_is_running_get(job_id=job_id)

        return GetJobStatusResponse(
            running=result.is_running
        )
