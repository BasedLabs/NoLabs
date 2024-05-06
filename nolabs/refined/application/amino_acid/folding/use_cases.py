__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature',
    'GetJobStatusFeature'
]

from typing import List, Tuple
from uuid import UUID

import esmfold_light_microservice
import esmfold_microservice
import rosettafold_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.amino_acid.folding.api_models import GetJobStatusResponse, JobResult, JobResponse, \
    SetupJobRequest
from nolabs.refined.domain.models.common import JobId, Experiment, JobName, Protein
from nolabs.refined.domain.models.folding import FoldingJob, FoldingBackendEnum
from nolabs.utils import generate_uuid


def map_job_to_response(job: FoldingJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        backend=FoldingBackendEnum(job.backend),
        proteins=[p.iid.value for p in job.proteins],
        result=[
            JobResult(
                protein_id=item.protein_id,
                pdb=item.pdb_content.decode('utf-8')
            )
            for item in job.foldings
        ]
    )


class GetJobFeature:
    """
    Use case - get job information.
    """
    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            job_id = JobId(job_id)
            job: FoldingJob = FoldingJob.objects.with_id(job_id.value)

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
            job_name = JobName(request.job_name if request.job_name else 'New folding job')
            folding_backend = FoldingBackendEnum(request.backend) if request.backend else FoldingBackendEnum.esmfold

            experiment = Experiment.objects.with_id(id=request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job: FoldingJob = FoldingJob.objects.with_id(job_id.value)

            if not job:
                job = FoldingJob(
                    id=job_id,
                    name=job_name,
                    experiment=experiment
                )

            proteins: List[Protein] = []
            for protein_id in request.proteins:
                protein = Protein.objects.get(id=protein_id)

                if not protein:
                    raise NoLabsException(ErrorCodes.protein_not_found)

                proteins.append(protein)

            job.set_inputs(proteins=proteins, backend=folding_backend)
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
    _esmfold: esmfold_microservice.DefaultApi
    _esmfold_light: esmfold_microservice.DefaultApi
    _rosettafold: rosettafold_microservice.DefaultApi

    def __init__(self,
                 esmfold: esmfold_microservice.DefaultApi,
                 esmfold_light: esmfold_microservice.DefaultApi,
                 rosettafold: rosettafold_microservice.DefaultApi):
        self._esmfold = esmfold
        self._esmfold_light = esmfold_light
        self._rosettafold = rosettafold

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        try:
            job: FoldingJob = FoldingJob.objects.with_id(job_id)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            folding_backend = job.backend

            if folding_backend == FoldingBackendEnum.esmfold:
                response = self._esmfold.is_job_running_job_job_id_is_running_get(job_id=job.iid.value)
                return GetJobStatusResponse(
                    running=response.is_running
                )

            if folding_backend == FoldingBackendEnum.esmfold_light:
                response = self._esmfold_light.is_job_running_job_job_id_is_running_get(job_id=job.iid.value)
                return GetJobStatusResponse(
                    running=response.is_running
                )

            if folding_backend == FoldingBackendEnum.rosettafold:
                response = self._rosettafold.is_job_running_job_job_id_is_running_get(job_id=job.iid.value)
                return GetJobStatusResponse(
                    running=response.is_running
                )

            raise NoLabsException(ErrorCodes.invalid_folding_backend)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e


class RunJobFeature:
    """
    Use case - start job.
    """
    _esmfold: esmfold_microservice.DefaultApi
    _esmfold_light: esmfold_microservice.DefaultApi
    _rosettafold: rosettafold_microservice.DefaultApi

    def __init__(self,
                 esmfold: esmfold_microservice.DefaultApi,
                 esmfold_light: esmfold_microservice.DefaultApi,
                 rosettafold: rosettafold_microservice.DefaultApi):
        self._esmfold = esmfold
        self._esmfold_light = esmfold_light
        self._rosettafold = rosettafold

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            assert job_id

            job_id = JobId(job_id)
            job: FoldingJob = FoldingJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            result: List[Tuple[Protein, str]] = []

            for protein in job.proteins:
                sequence = protein.get_fasta()
                pdb_content, errors = self.request_factory(job_id=job.iid, sequence=sequence,
                                                           folding_backend=FoldingBackendEnum(job.backend))

                if errors:
                    raise NoLabsException(ErrorCodes.amino_acid_localisation_run_error)

                result.append((protein, pdb_content))

            job.set_result(result=[
                (protein, pdb) for protein, pdb in result
            ])

            job.save(cascade=True)

            for protein, pdb in result:
                protein.set_pdb(pdb)
                protein.save()

            return map_job_to_response(job)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e

    def request_factory(self, job_id: JobId, sequence: str, folding_backend: FoldingBackendEnum) -> Tuple[
        str, List[str]]:
        if folding_backend == FoldingBackendEnum.esmfold:
            response = self._esmfold.predict_run_folding_post(
                run_esm_fold_prediction_request=esmfold_microservice.RunEsmFoldPredictionRequest(
                    protein_sequence=sequence
                ), _request_timeout=(1000.0, 1000.0))
            return (response.pdb_content, response.errors)

        if folding_backend == FoldingBackendEnum.esmfold_light:
            response = self._esmfold_light.predict_run_folding_post(
                run_esm_fold_prediction_request=esmfold_light_microservice.RunEsmFoldPredictionRequest(
                    protein_sequence=sequence
                ), _request_timeout=(1000.0, 1000.0))
            return (response.pdb_content, response.errors)

        if folding_backend == FoldingBackendEnum.rosettafold:
            response = self._rosettafold.run_folding_run_folding_post(
                job_id=job_id.value,
                fasta=sequence,
                _request_timeout=(1000.0, 1000.0)
            )
            return (response.pdb_content.anyof_schema_1_validator, response.errors)

        raise NoLabsException(ErrorCodes.folding_method_unknown, [str(folding_backend)])
