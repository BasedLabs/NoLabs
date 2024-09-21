from typing import List
from uuid import UUID

from domain.exceptions import ErrorCodes, NoLabsException

from nolabs.application.use_cases.jobs.api_models import (
    GetJobMetadataResponse, UpdateJobRequest)
from nolabs.domain.models.common import (Experiment, ExperimentId, Job, JobId,
                                         JobName)


class UpdateJobFeature:
    async def handle(self, job_id: UUID, request: UpdateJobRequest):
        job_id = JobId(job_id)
        job: Job = Job.objects.with_id(job_id.value)

        if not job:
            return

        if job.name != request.job_name:
            job.set_name(JobName(request.job_name))

        await job.save()


class DeleteJobFeature:
    async def handle(self, job_id: UUID):
        job_id = JobId(job_id)
        job: Job = Job.objects.with_id(job_id.value)
        if not job:
            return
        job.delete()


class GetJobsMetadataFeature:
    async def handle(self, experiment_id: UUID) -> List[GetJobMetadataResponse]:
        assert experiment_id

        experiment_id = ExperimentId(experiment_id)
        experiment = Experiment.objects.with_id(experiment_id.value)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        result: List[GetJobMetadataResponse] = []
        job: Job
        for job in Job.objects.filter(experiment=experiment):
            result.append(
                GetJobMetadataResponse(
                    job_id=job.id, job_name=str(job.name), type=str(type(job))
                )
            )

        return result


class GetJobMetadataFeature:
    async def handle(self, job_id: UUID) -> GetJobMetadataResponse:
        try:
            if not job_id:
                raise NoLabsException(ErrorCodes.job_not_found)

            job: Job = Job.objects.with_id(job_id)
            return GetJobMetadataResponse(
                job_id=job.id, job_name=str(job.name), type=str(type(job))
            )
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e

            raise NoLabsException(ErrorCodes.unknown_exception) from e
