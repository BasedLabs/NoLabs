__all__ = [
    'GetJobFeature',
    'SetupJobFeature',
    'RunJobFeature'
]

import asyncio
import uuid
from typing import List
from uuid import UUID

from mongoengine import Q

from nolabs.application.rfdiffusion.api_models import JobResponse, SetupJobRequest
from nolabs.application.rfdiffusion.worker_models import RunRfdiffusionRequest, RunRfdiffusionResponse
from nolabs.application.rfdiffusion.workflow import RfDiffusionInput, RfDiffusionOutput
from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import JobId, JobName, Experiment, Protein, ProteinName
from nolabs.domain.models.protein_design import RfdiffusionJob
from nolabs.infrastructure.celery_app_factory import get_celery_app, wait_for_task


def map_job_to_response(job: RfdiffusionJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        protein_id=job.protein.iid.value,
        binder_ids=[p.iid.value for p in job.binders],
        number_of_designs=job.number_of_designs,
        timesteps=job.timesteps,
        hotspots=job.hotspots,
        contig=job.contig
    )


class GetJobFeature:
    """
    Use case - get job information.
    """
    async def handle(self, job_id: UUID) -> JobResponse:
        job_id = JobId(job_id)
        job: RfdiffusionJob = RfdiffusionJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return map_job_to_response(job)


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """
    async def handle(self, request: SetupJobRequest) -> JobResponse:
        assert request

        job_id = JobId(request.job_id if request.job_id else uuid.uuid4())
        job_name = JobName(request.job_name if request.job_name else 'New protein design job')

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job: RfdiffusionJob = RfdiffusionJob.objects(Q(id=job_id.value) | Q(name=job_name.value)).first()

        if not job:
            job = RfdiffusionJob(
                id=job_id,
                name=job_name,
                experiment=experiment
            )

        protein = Protein.objects.with_id(request.protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        job.set_input(
            protein=protein,
            contig=request.contig,
            number_of_designs=request.number_of_designs,
            hotspots=request.hotspots,
            timesteps=request.timesteps
        )
        await job.save(cascade=True)

        return map_job_to_response(job)


class RunJobFeature:
    """
    Use case - start job.
    """
    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: RfdiffusionJob = RfdiffusionJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        task_id = get_celery_app().send_task(name='design', queue='rfdiffusion',
                                             kwargs={'param': RunRfdiffusionRequest(
                                                pdb_content=job.protein.get_pdb(),
                                                 contig=job.contig,
                                                 hotspots=job.hotspots,
                                                 timesteps=job.timesteps,
                                                 number_of_designs=job.number_of_designs
                                             )})

        job_result = await wait_for_task(task_id=task_id)
        job_result = RunRfdiffusionResponse(**job_result)

        if job_result.errors and not job_result.pdbs_content:
            raise NoLabsException(ErrorCodes.protein_design_run_error, ', '.join(job_result.errors))

        binders: List[Protein] = []
        for i, pdb in enumerate(job_result.pdbs_content):
            binder = Protein.create(
                experiment=job.component.experiment.id,
                name=ProteinName(f'{job.protein.name.value}-binder'),
                pdb_content=pdb
            )
            binders.append(binder)

            job.protein.add_protein_binder(binder)
            await job.save(cascade=True)
            binder.save(cascade=True)

        job.set_result(binders=binders)

        return map_job_to_response(job)

    async def design(
        self, task_id: str, payload: RfDiffusionInput
    ) -> RfDiffusionOutput:
        app = get_celery_app()
        async_result = app.send_task(
            id=task_id,
            name="design",
            queue="rfdiffusion",
            args=[payload.model_dump()],
        )

        while not async_result.ready():
            await asyncio.sleep(0.5)
        result = async_result.get()
        return RfDiffusionOutput(**result)
