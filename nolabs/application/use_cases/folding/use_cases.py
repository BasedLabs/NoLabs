__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature',
    'GetJobStatusFeature'
]

import tempfile
from typing import List, Tuple, Optional
from uuid import UUID

import esmfold_light_microservice
import esmfold_microservice
import rosettafold_microservice
from mongoengine import Q

from domain.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.folding.api_models import GetJobStatusResponse, JobResult, JobResponse, \
    SetupJobRequest
from nolabs.domain.models.common import JobId, Experiment, JobName, Protein
from nolabs.domain.models.folding import FoldingJob, FoldingBackendEnum
from nolabs.utils import generate_uuid


def map_job_to_response(job: FoldingJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        backend=FoldingBackendEnum(job.backend),
        protein_ids=[job.protein.id],
        result=[
            JobResult(
                protein_id=job.folding.protein_id,
                pdb=job.folding.pdb_content.decode('utf-8')
            )
        ] if job.folding else [],
        experiment_id=job.experiment.id
    )


class GetJobFeature:
    """
    Use case - get job information.
    """
    async def handle(self, job_id: UUID) -> JobResponse:
        job_id = JobId(job_id)
        job: FoldingJob = FoldingJob.objects.with_id(job_id.value)

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
        job_name = JobName(request.job_name if request.job_name else 'New folding job')
        folding_backend = FoldingBackendEnum(request.backend) if request.backend else FoldingBackendEnum.esmfold

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job: FoldingJob = FoldingJob.objects(Q(id=job_id.value) | Q(name=job_name.value)).first()

        if not job:
            job = FoldingJob(
                id=job_id,
                name=job_name,
                experiment=experiment
            )

        proteins: List[Protein] = []
        for protein_id in request.protein_ids:
            protein = Protein.objects.with_id(protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            proteins.append(protein)

        job.set_inputs(proteins=proteins, backend=folding_backend)
        await job.save(cascade=True)

        return map_job_to_response(job)


class GetJobStatusFeature:
    """
    Use case - set job status.
    """
    _esmfold: esmfold_microservice.DefaultApi
    _esmfold_light: esmfold_light_microservice.DefaultApi
    _rosettafold: rosettafold_microservice.DefaultApi

    def __init__(self,
                 esmfold: esmfold_microservice.DefaultApi,
                 esmfold_light: esmfold_light_microservice.DefaultApi,
                 rosettafold: rosettafold_microservice.DefaultApi):
        self._esmfold = esmfold
        self._esmfold_light = esmfold_light
        self._rosettafold = rosettafold

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        return GetJobStatusResponse(running=True, result_valid=True) # TODO remove

        job: FoldingJob = FoldingJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        folding_backend = job.backend

        if folding_backend == FoldingBackendEnum.esmfold:
            response = self._esmfold.is_job_running_job_job_id_is_running_get(job_id=str(job.iid.value))
            return GetJobStatusResponse(
                running=response['is_running'],
                result_valid=job.result_valid()
            )

        if folding_backend == FoldingBackendEnum.esmfold_light:
            response = self._esmfold_light.is_job_running_job_job_id_is_running_get(job_id=str(job.iid.value))
            return GetJobStatusResponse(
                running=response.is_running,
                result_valid=job.result_valid()
            )

        if folding_backend == FoldingBackendEnum.rosettafold:
            response = self._rosettafold.is_job_running_job_job_id_is_running_get(job_id=str(job.iid.value))
            return GetJobStatusResponse(
                running=response.is_running,
                result_valid=job.result_valid()
            )

        raise NoLabsException(ErrorCodes.invalid_folding_backend)


class RunJobFeature:
    """
    Use case - start job.
    """
    _esmfold: esmfold_microservice.DefaultApi
    _esmfold_light: esmfold_microservice.DefaultApi
    _rosettafold: rosettafold_microservice.DefaultApi

    def __init__(self,
                 esmfold: Optional[esmfold_microservice.DefaultApi] = None,
                 esmfold_light: Optional[esmfold_light_microservice.DefaultApi] = None,
                 rosettafold: Optional[rosettafold_microservice.DefaultApi] = None):
        if not esmfold and not esmfold_light and not rosettafold:
            raise NoLabsException(ErrorCodes.folding_run_error, 'You must specify at least one folding backend')

        self._esmfold = esmfold
        self._esmfold_light = esmfold_light
        self._rosettafold = rosettafold

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: FoldingJob = FoldingJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        result: List[Tuple[Protein, str]] = []
        try:
            job.started()
            await job.save()


            for protein in job.proteins:
                sequence = protein.get_amino_acid_sequence()

                if not sequence:
                    raise NoLabsException(ErrorCodes.protein_fasta_is_empty, f'Protein sequence is empty for protein {protein.name}')

                pdb_content, errors = self.request_factory(job_id=job.iid, sequence=sequence,
                                                           folding_backend=FoldingBackendEnum(job.backend))

                if errors:
                    raise NoLabsException(ErrorCodes.folding_run_error, errors)

                result.append((protein, pdb_content))

            job.set_result(result=[
                (protein, pdb) for protein, pdb in result
            ])

            await job.save(cascade=True)

            for protein, pdb in result:
                protein.set_pdb(pdb)
                protein.save()

        finally:
            job.finished()
            await job.save()

        return map_job_to_response(job)

    def request_factory(self, job_id: JobId, sequence: str, folding_backend: FoldingBackendEnum) -> Tuple[
        str, List[str]]:
        if folding_backend == FoldingBackendEnum.esmfold:
            response = self._esmfold.predict_run_folding_post(
                run_esm_fold_prediction_request=esmfold_microservice.RunEsmFoldPredictionRequest(
                    protein_sequence=sequence
                ), _request_timeout=(1000.0, 1000.0))
            return (response.pdb_content, response.errors)

        if folding_backend == FoldingBackendEnum.esmfold_light:
            response = self._esmfold_light.predict_through_api_run_folding_post(
                run_esm_fold_prediction_request=esmfold_light_microservice.RunEsmFoldPredictionRequest(
                    protein_sequence=sequence
                ), _request_timeout=(1000.0, 1000.0))
            return (response.pdb_content, [])

        if folding_backend == FoldingBackendEnum.rosettafold:
            with tempfile.NamedTemporaryFile(mode="w+t", delete=False) as temp_file:
                temp_file.write(sequence)
                temp_file.close()
                response = self._rosettafold.run_folding_run_folding_post(
                    fasta=temp_file.name,
                    _request_timeout=(1000.0, 1000.0)
                )
                return (response.pdb_content.anyof_schema_1_validator, response.errors)

        raise NoLabsException(ErrorCodes.folding_method_unknown, [str(folding_backend)])
