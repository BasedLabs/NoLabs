__all__ = [
    'GetJobStatusFeature',
    'RunMsaGenerationJobFeature',
    'GetMsaPredictionJobResultFeature'
]


from uuid import UUID

import msa_light_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.msa_generation.api_models import PredictMsaRequest, PredictMsaResponse, \
    GetJobStatusResponse
from nolabs.refined.domain.models import Experiment, JobId, JobName, Protein, Job
from nolabs.refined.domain.models.msa import MsaGenerationJob
from nolabs.refined.infrastructure.settings import MsaLightMicroserviceSettings


class RunMsaGenerationJobFeature:
    def __init__(self, settings: MsaLightMicroserviceSettings, api: msa_light_microservice.DefaultApi):
        self._api = api
        self._settings = settings

    async def handle(self, request: PredictMsaRequest) -> PredictMsaResponse:
        assert request

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job = MsaGenerationJob.objects.with_id(request.job_id)

        if not job:
            job = MsaGenerationJob(
                id=JobId(request.job_id),
                name=JobName('New msa generation job'),
                experiment=experiment)

        protein = Protein.objects.with_id(request.protein_id)

        if not protein.get_fasta():
            raise NoLabsException(ErrorCodes.protein_fasta_is_empty,
                                  f'No fasta content for the protein with id {protein.iid.value}')

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        job.set_input(protein=protein)

        request = msa_light_microservice.RunMsaPredictionRequest(
            api_url=self._settings.msa_server_url,
            fasta_contents=protein.get_fasta())
        msa_contents = self._api.predict_msa_predict_msa_post(predict_msa_predict_msa_post=request).msa_contents

        job.set_result(protein=protein, msa=msa_contents)
        protein.set_msa(msa=msa_contents)

        return PredictMsaResponse(msa_contents)


class GetJobStatusFeature:
    def __init__(self, api: msa_light_microservice.DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        job_id = JobId(job_id)

        if not Job.objects.with_id(job_id):
            raise NoLabsException(ErrorCodes.job_not_found)

        response = self._api.is_job_running_job_job_id_is_running_get(job_id=job_id)

        return GetJobStatusResponse(running=response.is_running)


class GetMsaPredictionJobResultFeature:
    async def handle(self, job_id: UUID) -> PredictMsaResponse:
        job_id = JobId(job_id)

        job: MsaGenerationJob = MsaGenerationJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return PredictMsaResponse(
            msa=job.msa
        )
