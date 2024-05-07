__all__ = [
    'GetJobStatusFeature',
    'RunJobFeature',
    'GetJobFeature',
    'SetupJobFeature'
]

from uuid import UUID

import msa_light_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.msa_generation.api_models import JobResponse, SetupJobRequest, GetJobStatusResponse
from nolabs.refined.domain.models import Experiment, JobId, JobName, Protein, Job
from nolabs.refined.domain.models.msa import MsaGenerationJob
from nolabs.refined.infrastructure.settings import MsaLightMicroserviceSettings
from nolabs.utils import generate_uuid


def map_job_to_response(job: MsaGenerationJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        experiment_id=job.experiment.iid.value,
        protein_id=job.protein.iid.value,
        msa=job.msa.decode('utf-8') if job.msa else None
    )


class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            job_id = JobId(job_id)
            job: MsaGenerationJob = MsaGenerationJob.objects.with_id(job_id.value)

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
            job_name = JobName(request.job_name if request.job_name else 'New msa job')

            experiment = Experiment.objects.with_id(id=request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job: MsaGenerationJob = MsaGenerationJob.objects.with_id(job_id.value)

            if not job:
                job = MsaGenerationJob(
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


class RunJobFeature:
    """
    Use case - start job.
    """

    def __init__(self, settings: MsaLightMicroserviceSettings, api: msa_light_microservice.DefaultApi):
        self._api = api
        self._settings = settings

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            assert job_id

            job = MsaGenerationJob.objects.with_id(job_id)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            request = msa_light_microservice.RunMsaPredictionRequest(
                api_url=self._settings.msa_server_url,
                fasta_contents=job.protein.get_fasta())
            msa_contents = self._api.predict_msa_predict_msa_post(predict_msa_predict_msa_post=request).msa_contents

            job.set_result(protein=job.protein, msa=msa_contents)
            job.save(cascade=True)
            job.protein.set_msa(msa=msa_contents)
            job.protein.save()

            return map_job_to_response(job)
        except Exception as e:
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e


class GetJobStatusFeature:
    def __init__(self, api: msa_light_microservice.DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        try:
            job_id = JobId(job_id)

            job: MsaGenerationJob = MsaGenerationJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            response = self._api.is_job_running_job_job_id_is_running_get(job_id=job_id)

            return GetJobStatusResponse(
                running=response.is_running
            )
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e
