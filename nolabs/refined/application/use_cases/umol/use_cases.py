__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature',
    'GetJobStatusFeature'
]

from uuid import UUID

import umol_microservice
from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.umol.api_models import JobResponse, JobResult, \
    SetupJobRequest, GetJobStatusResponse
from nolabs.refined.domain.models.common import JobId, Experiment, JobName, Protein, Ligand
from nolabs.refined.domain.models.umol import UmolBindingJob
from nolabs.utils import generate_uuid


def map_job_to_response(job: UmolBindingJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        pocket_ids=job.pocket_ids,
        protein_id=job.protein.iid.value,
        ligand_id=job.ligand.iid.value,
        result=JobResult(
            predicted_pdb=job.predicted_pdb,
            predicted_sdf=job.predicted_sdf,
            plddt_array=job.plddt_array
        )
    )


class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            job_id = JobId(job_id)
            job: UmolBindingJob = UmolBindingJob.objects.with_id(job_id.value)

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
            job_name = JobName(request.job_name if request.job_name else 'New protein ligand UMOL binding job')

            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job: UmolBindingJob = UmolBindingJob.objects(Q(id=job_id.value) | Q(name=job_name.value)).first()

            if not job:
                job = UmolBindingJob(
                    id=job_id,
                    name=job_name,
                    experiment=experiment
                )

            protein = Protein.objects(id=request.protein_id, experiment=experiment).first()

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            ligand = Ligand.objects(id=request.protein_id, experiment=experiment).first()

            job.set_input(
                protein=protein,
                ligand=ligand,
                pocket_ids=request.pocket_ids
            )

            job.save(cascade=True)

            return map_job_to_response(job)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e


class GetJobStatusFeature:
    """
    Use case - set job status.
    """
    _umol = umol_microservice.DefaultApi

    def __init__(self,
                 umol: umol_microservice.DefaultApi):
        self._umol = umol

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        try:
            job: UmolBindingJob = UmolBindingJob.objects.with_id(job_id)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            response = self._umol.is_job_running_job_job_id_is_running_get(
                job_id=job.iid.value
            )
            return GetJobStatusResponse(
                running=response.is_running
            )
        except Exception as e:
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e


class RunJobFeature:
    """
    Use case - start job.
    """
    _umol = umol_microservice.DefaultApi

    def __init__(self, umol: umol_microservice.DefaultApi):
        self._umol = umol

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            assert job_id

            job_id = JobId(job_id)
            job: UmolBindingJob = UmolBindingJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            request = umol_microservice.RunUmolPredictionRequest(
                protein_sequence=job.protein.get_amino_acid_sequence(),
                ligand_smiles=job.ligand.get_smiles(),
                pocket_ids=job.pocket_ids,
                job_id=job_id.value,
                msa_content=job.protein.get_msa()
            )
            response = self._umol.predict_run_umol_post(run_umol_prediction_request=request)

            if response.errors:
                raise NoLabsException(ErrorCodes.diffdock_api_error, response.errors)

            predicted_pdb = response.pdb_contents
            predicted_sdf = response.sdf_contents
            plddt_array = response.plddt_array

            job.set_result(
                protein=job.protein,
                ligand=job.ligand,
                predicted_sdf=predicted_sdf,
                predicted_pdb=predicted_pdb,
                plddt_array=plddt_array
            )
            job.save(cascade=True)

            ligand = job.ligand
            ligand.add_binding(
                protein=job.protein,
                sdf_content=predicted_sdf,
                minimized_affinity=None,
                scored_affinity=None,
                confidence=None,
                plddt_array=plddt_array,
                pdb_content=predicted_pdb
            )
            ligand.save(cascade=True)

            return map_job_to_response(job)

        except Exception as e:
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e
