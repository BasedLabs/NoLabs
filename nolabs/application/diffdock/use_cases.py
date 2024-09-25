__all__ = ["GetJobFeature", "RunJobFeature", "SetupJobFeature"]

from typing import List
from uuid import UUID, uuid4

from microservices.diffdock.service.api_models import RunDiffDockPredictionRequest

from nolabs.application.diffdock.api_models import (JobResponse, JobResult, SetupJobRequest)
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import (Experiment, JobId, JobName, Ligand,
                                         Protein)
from nolabs.domain.models.diffdock import DiffDockBindingJob
from nolabs.infrastructure.cel import cel as celery
from nolabs.utils import generate_uuid


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

        job_id = JobId(request.job_id if request.job_id else generate_uuid())
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

        response = celery.diffdock_inference(
            task_id=task_id,
            payload=RunDiffDockPredictionRequest(
                pdb_contents=job.protein.get_pdb(),
                sdf_contents=ligand.get_sdf(),
                samples_per_complex=job.samples_per_complex,
            ),
        )

        for item in response.sdf_results:
            complex = ligand.add_binding(
                protein=job.protein,
                sdf_content=item.sdf_content,
                minimized_affinity=item.minimized_affinity,
                scored_affinity=item.scored_affinity,
                confidence=item.confidence,
                plddt_array=[],
                name=item.sdf_file_name,
            )
            result.append(complex)

        job.set_result(complexes=result)
        await job.save(cascade=True)

        return map_job_to_response(job)
