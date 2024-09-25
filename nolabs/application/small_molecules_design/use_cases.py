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

import uuid
from typing import List
from uuid import UUID

from domain.exceptions import ErrorCodes, NoLabsException
from microservices.reinvent.service.api_models import (
    RunReinforcementLearningRequest, RunSamplingRequest)

from nolabs.application.small_molecules_design.api_models import (
    GetJobStatusResponse, JobResponse, LogsResponse, SetupJobRequest,
    SmilesResponse, StartSamplingRequest)
from nolabs.application.small_molecules_design.services import \
    ReinventParametersSaver
from nolabs.domain.models.common import (Experiment, Job, JobId, JobName,
                                         Protein)
from nolabs.domain.models.small_molecules_design import SmallMoleculesDesignJob
from nolabs.infrastructure.cel import cel as celery
from nolabs.infrastructure.settings import settings
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


reinvent_directory = settings.reinvent_directory


class DeleteJobFeature:
    async def handle(self, job_id: UUID):
        (reinvent_directory / str(job_id)).rmdir()
        job_id = JobId(job_id)
        job: Job = Job.objects.with_id(job_id.value)
        if not job:
            return
        job.delete()


class GetJobStatusFeature:
    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(job_id)

        if not job:
            return GetJobStatusResponse(
                running=False, sampling_allowed=False, result_valid=False
            )

        chkpt_exists = (reinvent_directory / str(job_id) / "rl_direct.chkpt").exists()

        celery_task_ready = False
        celery_task_id = job.celery_task_id

        if celery_task_id:
            task_result = celery.task_result(celery_task_id)
            celery_task_ready = task_result.ready()

        return GetJobStatusResponse(
            running=celery_task_ready,
            sampling_allowed=not chkpt_exists and not celery_task_ready,
            result_valid=job.result_valid(),
        )


class GetJobFeature:
    async def handle(self, job_id: UUID) -> JobResponse:
        if not Job.objects(id=job_id):
            raise NoLabsException(ErrorCodes.job_not_found)

        job_id = JobId(job_id)
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(
            job_id.value
        )

        return map_job_to_response(job)


class GetJobLogsFeature:
    async def handle(self, job_id: UUID) -> LogsResponse:
        if not Job.objects(id=job_id):
            return LogsResponse(output="", docking_output="", errors="")

        output = (reinvent_directory / str(job_id) / "output.log").read_text()
        errors = (reinvent_directory / str(job_id) / "error.log").read_text()
        docking_output = (reinvent_directory / str(job_id) / "docking.log").read_text()

        return LogsResponse(
            output=output,
            docking_output=docking_output,
            errors=errors,
        )


class GetJobSmilesFeature:
    async def handle(self, job_id: UUID) -> List[SmilesResponse]:
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(job_id)
        if not job:
            return []

        ligands = job.ligands

        result = []

        for ligand in ligands:
            result.append(
                SmilesResponse(
                    smiles=ligand.get_smiles(),
                    drug_likeness=ligand.drug_likeness.value,
                    score=ligand.designed_ligand_score.value,
                    created_at=ligand.created_at,
                    stage=ligand.generated_stage,
                )
            )

        return result


class SetupJobFeature:
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
                message="Protein pdb content is undefined",
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

        parameters_saver = ReinventParametersSaver()
        await parameters_saver.save_params(job=job, pdb=job.protein.pdb_content)

        job.change_sampling_size(request.sampling_size)

        await job.save(cascade=True)

        return map_job_to_response(job)


class RunLearningStageJobFeature:
    async def handle(self, job_id: UUID):
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        if not job.input_errors():
            task_id = uuid.uuid4()
            await celery.reinvent_run_learning(
                task_id=task_id,
                request=RunReinforcementLearningRequest(config_id=str(job_id)),
            )
            job.set_task_id(task_id=task_id)
            await job.save()
            return

        input_error = job.input_errors()[0]
        raise NoLabsException(
            error_code=input_error.error_code, message=input_error.message
        )


class RunSamplingStageJobFeature:
    async def handle(self, job_id: UUID, request: StartSamplingRequest):
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        chkpt_exists = (reinvent_directory / str(job_id) / "rl_direct.chkpt").exists()

        celery_task_ready = False
        celery_task_id = job.celery_task_id

        if celery_task_id:
            task_result = celery.task_result(celery_task_id)
            celery_task_ready = task_result.ready()

        sampling_allowed = not chkpt_exists and not celery_task_ready

        if not sampling_allowed:
            raise NoLabsException(
                ErrorCodes.reinvent_cannot_run_sampling,
                "Cannot run sampling because learning stage was not started first or sampling is already running",
            )

        (res, task_id) = await celery.reinvent_run_sampling(
            request=RunSamplingRequest(
                config_id=str(job_id),
                number_of_molecules_to_generate=request.sampling_size,
            ),
            wait=False,
        )
        job.set_task_id(task_id=task_id)
        await job.save()


class StopJobFeature:
    async def handle(self, job_id: UUID):
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        if job.celery_task_id:
            celery.cancel_task(job.celery_task_id)
