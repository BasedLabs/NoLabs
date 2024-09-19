__all__ = [
    "DeleteJobFeature",
    "GetJobStatusFeature",
    "GetJobFeature",
    "GetJobLogsFeature",
    "GetJobSmilesFeature",
    "SetupJobFeature",
    "RunLearningStageJobFeature",
    "RunSamplingStageJobFeature",
    "StopJobFeature",
]

import os
from typing import List
from uuid import UUID

import reinvent_microservice
from domain.exceptions import ErrorCodes, NoLabsException

from nolabs.application.use_cases.small_molecules_design.api_models import (
    GetJobStatusResponse, JobResponse, LogsResponse, SetupJobRequest,
    SmilesResponse)
from nolabs.domain.models.common import (Experiment, Job, JobId, JobName,
                                         Protein)
from nolabs.domain.models.small_molecules_design import SmallMoleculesDesignJob
from nolabs.utils import generate_uuid


def map_job_to_response(job: SmallMoleculesDesignJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        experiment_id=job.experiment.iid.value,
        protein_id=job.protein.iid.value,
        center_x=job.center_x,
        center_y=job.center_y,
        center_z=job.center_z,
        size_x=job.size_x,
        size_y=job.size_y,
        size_z=job.size_z,
        batch_size=job.batch_size,
        minscore=job.minscore,
        epochs=job.epochs,
        sampling_size=job.sampling_size,
    )


class DeleteJobFeature:
    async def handle(self, job_id: UUID):
        job_id = JobId(job_id)
        job: Job = Job.objects.with_id(job_id.value)
        if not job:
            return
        job.delete()


class GetJobStatusFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        if not SmallMoleculesDesignJob.objects(id=job_id):
            return GetJobStatusResponse(
                running=False, sampling_allowed=False, result_valid=False
            )

        job = SmallMoleculesDesignJob.objects.with_id(job_id)

        config_result = self._api.get_config_api_reinvent_reinvent_config_id_get(
            config_id=str(job_id)
        )

        if not config_result:
            return GetJobStatusResponse(
                running=False, sampling_allowed=False, result_valid=job.result_valid()
            )

        config = config_result.actual_instance

        return GetJobStatusResponse(
            running=config.running,
            sampling_allowed=config.sampling_allowed,
            result_valid=job.result_valid(),
        )


class GetJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID) -> JobResponse:
        if not Job.objects(id=job_id):
            raise NoLabsException(ErrorCodes.job_not_found)

        job_id = JobId(job_id)
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(
            job_id.value
        )

        return map_job_to_response(job)


class GetJobLogsFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID) -> LogsResponse:
        if not Job.objects(id=job_id):
            return LogsResponse(output="", docking_output="", errors="")

        response = self._api.logs_api_reinvent_config_id_logs_get(config_id=str(job_id))

        if not response:
            return LogsResponse(output="Empty", docking_output="Empty", errors="Empty")

        instance = response.actual_instance

        return LogsResponse(
            output=instance.output,
            docking_output=instance.docking_output,
            errors=instance.errors,
        )


class GetJobSmilesFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID) -> List[SmilesResponse]:
        job = Job.objects.with_id(job_id)
        if not job:
            return []

        response = self._api.smiles_api_reinvent_config_id_smiles_get(
            config_id=str(job_id)
        )

        return [
            SmilesResponse(
                smiles=s.smiles,
                drug_likeness=s.drug_likeness,
                score=s.score,
                stage=s.stage,
                created_at=job.created_at,
            )
            for s in response.smiles
        ]


class SetupJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, request: SetupJobRequest) -> JobResponse:
        job_id = JobId(request.job_id if request.job_id else generate_uuid())
        job_name = JobName(
            request.job_name if request.job_name else "New small molecules design job"
        )

        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(
            job_id.value
        )

        if not job:
            if not request.experiment_id:
                raise NoLabsException(ErrorCodes.invalid_experiment_id)

            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job = SmallMoleculesDesignJob(
                id=JobId(request.job_id), name=job_name, experiment=experiment
            )

        protein = Protein.objects.with_id(request.protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        if not protein.pdb_content:
            raise NoLabsException(
                ErrorCodes.protein_pdb_is_empty,
                messages="Protein pdb content is undefined",
            )

        job.set_inputs(
            protein=protein,
            center_x=request.center_x,
            center_y=request.center_y,
            center_z=request.center_z,
            size_x=request.size_x,
            size_y=request.size_y,
            size_z=request.size_z,
            batch_size=request.batch_size,
            minscore=request.minscore,
            epochs=request.epochs,
        )

        tmp_file_path = "tmp.pdb"
        open(tmp_file_path, "wb").write(job.protein.pdb_content)

        self._api.save_params_api_reinvent_config_id_params_post(
            config_id=str(job.iid.value),
            name=job.name.value,
            pdb_file=os.path.abspath(tmp_file_path),
            center_x=request.center_x,
            center_y=request.center_y,
            center_z=request.center_z,
            size_x=request.size_x,
            size_y=request.size_y,
            size_z=request.size_z,
            batch_size=request.batch_size,
            minscore=request.minscore,
            epochs=request.epochs,
        )

        job.change_sampling_size(request.sampling_size)

        await job.save(cascade=True)

        return map_job_to_response(job)


class RunLearningStageJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID):
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        if not job.input_errors():
            self._api.learning_api_reinvent_config_id_start_learning_post(
                config_id=str(job_id)
            )
            return

        input_error = job.input_errors()[0]
        raise NoLabsException(
            error_code=input_error.error_code, messages=input_error.message
        )


class RunSamplingStageJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID):
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        config_result = self._api.get_config_api_reinvent_reinvent_config_id_get(
            config_id=str(job_id)
        )

        if not config_result.actual_instance.sampling_allowed:
            raise NoLabsException(
                ErrorCodes.reinvent_cannot_run_sampling,
                "Cannot run sampling because learning stage was not started first or sampling is already running",
            )

        self._api.sampling_api_reinvent_config_id_start_sampling_post(
            config_id=str(job_id),
            sampling_size_request=reinvent_microservice.SamplingSizeRequest(
                number_of_molecules_to_design=job.sampling_size
            ),
        )


class StopJobFeature:
    def __init__(self, api: reinvent_microservice.ReinventApi):
        self._api = api

    async def handle(self, job_id: UUID):
        self._api.stop_api_reinvent_config_id_jobs_stop_post(config_id=str(job_id))
