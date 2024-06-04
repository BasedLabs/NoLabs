__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature',
    'GetJobStatusFeature'
]

from typing import List
from uuid import UUID

import diffdock_microservice
from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.diffdock.api_models import JobResponse, JobResult, \
    SetupJobRequest, GetJobStatusResponse
from nolabs.refined.domain.models.common import JobId, Experiment, JobName, Protein, Ligand
from nolabs.refined.domain.models.diffdock import DiffDockBindingJob, DiffDockJobResult
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
                complex_id=res.complex_id,
                sdf_content=res.sdf_content.decode('utf-8'),
                minimized_affinity=res.minimized_affinity,
                scored_affinity=res.scored_affinity,
                confidence=res.confidence
            )
            for res in job.result
        ]
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
        job_name = JobName(request.job_name if request.job_name else 'New protein ligand DIFFDOCK binding job')

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job: DiffDockBindingJob = DiffDockBindingJob.objects(Q(id=job_id.value) | Q(name=job_name.value)).first()

        if not job:
            job = DiffDockBindingJob(
                id=job_id,
                name=job_name,
                experiment=experiment
            )

        protein = Protein.objects(id=request.protein_id, experiment=experiment).first()

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        ligand = Ligand.objects(id=request.ligand_id, experiment=experiment).first()
        if not ligand:
            raise NoLabsException(ErrorCodes.ligand_is_undefined)

        job.set_input(protein=protein,
                      ligand=ligand,
                      samples_per_complex=request.samples_per_complex)
        job.save(cascade=True)

        return map_job_to_response(job)


class GetJobStatusFeature:
    """
    Use case - set job status.
    """
    _diffdock = diffdock_microservice.DefaultApi

    def __init__(self,
                 diffdock: diffdock_microservice.DefaultApi):
        self._diffdock = diffdock

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        response = self._diffdock.is_job_running_job_job_id_is_running_get(
            job_id=job.iid.value
        )
        return GetJobStatusResponse(
            running=response.is_running
        )


class RunJobFeature:
    """
    Use case - start job.
    """
    _diffdock = diffdock_microservice.DefaultApi

    def __init__(self, diffdock: diffdock_microservice.DefaultApi):
        self._diffdock = diffdock

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        result: List[DiffDockJobResult] = []

        ligand = job.ligand

        request = diffdock_microservice.RunDiffDockPredictionRequest(pdb_contents=job.protein.get_pdb(),
                                                                     sdf_contents=ligand.get_sdf(),
                                                                     samples_per_complex=job.samples_per_complex,
                                                                     job_id=str(job_id.value),
                                                                     )
        response = self._diffdock.predict_run_docking_post(
            run_diff_dock_prediction_request=request, _request_timeout=(1000.0, 1000.0))

        if not response.success:
            raise NoLabsException(ErrorCodes.diffdock_api_error, response.message)

        for item in response.sdf_results:
            complex = ligand.add_binding(
                protein=job.protein,
                sdf_content=item.sdf_content,
                minimized_affinity=item.minimized_affinity,
                scored_affinity=item.scored_affinity,
                confidence=item.confidence,
                plddt_array=[]
            )
            result.append(
                DiffDockJobResult(
                    complex_id=complex.iid.value,
                    sdf_content=item.sdf_content,
                    minimized_affinity=item.minimized_affinity,
                    scored_affinity=item.scored_affinity,
                    confidence=item.confidence
                )
            )

        job.set_result(protein=job.protein,
                       ligand=ligand,
                       result=result)
        job.save(cascade=True)

        return map_job_to_response(job)
