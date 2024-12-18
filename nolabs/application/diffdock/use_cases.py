__all__ = ["GetJobFeature", "RunJobFeature", "SetupJobFeature"]

from typing import List
from uuid import UUID, uuid4

from nolabs.application.diffdock.worker_models import RunDiffDockPredictionInput, RunDiffDockPredictionOutput
from nolabs.infrastructure.celery_app_factory import get_celery_app, wait_for_task
from nolabs.infrastructure.mongo_connector import get_connection

from nolabs.application.diffdock.api_models import (
    JobResponse,
    JobResult,
    SetupJobRequest,
)
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Experiment, JobId, JobName, Ligand, Protein
from nolabs.domain.models.diffdock import DiffDockBindingJob


def map_job_to_response(job: DiffDockBindingJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        samples_per_complex=job.samples_per_complex,
        protein_id=job.protein.iid.value,
        ligand_id=job.ligand.iid.value,
        result=[
            JobResult(
                complex_id=res.id,
                sdf_content=res.binding_ligand.get_sdf(),
                minimized_affinity=res.minimized_affinity,
                scored_affinity=res.scored_affinity,
                confidence=res.confidence,
            )
            for res in job.complexes
        ],
    )


class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        job_id = JobId(job_id)
        job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id.value)

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
        job_name = JobName(
            request.job_name
            if request.job_name
            else "New protein ligand DIFFDOCK binding job"
        )

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id.value)

        if not job:
            job = DiffDockBindingJob(id=job_id, name=job_name, experiment=experiment)

        protein = Protein.objects.with_id(request.protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        ligand = Ligand.objects(id=request.ligand_id, experiment=experiment).first()
        if not ligand:
            raise NoLabsException(ErrorCodes.ligand_is_undefined)

        job.set_input(
            protein=protein,
            ligand=ligand,
            samples_per_complex=request.samples_per_complex,
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
        job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        result: List[Protein] = []

        ligand = job.ligand

        task_id = uuid4()
        job.set_task_id(task_id=str(task_id))
        await job.save()

        task_id = get_celery_app().send_task(name='inference', queue='diffdock', kwargs={'param': RunDiffDockPredictionInput(
            pdb_contents=job.protein.get_pdb(),
            sdf_contents=job.ligand.get_sdf()
        )})

        job_result = await wait_for_task(task_id=task_id)
        job_result = RunDiffDockPredictionOutput(**job_result)

        db = get_connection()
        session = db.client.start_session()
        with session.start_transaction():
            for item in job_result.sdf_results:
                ligand_for_complex = Ligand.copy(ligand)

                ligand_for_complex.save(session=session)

                complex = Protein.create_complex(
                    protein=job.protein,
                    ligand=ligand_for_complex,
                    minimized_affinity=item.minimized_affinity,
                    scored_affinity=item.scored_affinity,
                    confidence=item.confidence,
                    plddt_array=[]
                )

                complex.save(session=session)
                result.append(complex)

            job.set_result(complexes=result)

            await job.save(session=session, cascade=True)

        job.set_result(complexes=result)
        await job.save(cascade=True)

        return map_job_to_response(job)