__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature'
]

from typing import List, Tuple
from uuid import UUID

from mongoengine import Q
from solubility_microservice import DefaultApi, RunSolubilityPredictionRequest

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.solubility.api_models import JobResponse, JobResult, SetupJobRequest, \
    GetJobStatusResponse
from nolabs.domain.models.common import JobId, Experiment, JobName, Protein, SolubleProbability
from nolabs.domain.models.solubility import SolubilityJob
from nolabs.utils import generate_uuid


def map_job_to_response(job: SolubilityJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        protein_ids=[p.iid.value for p in job.proteins],
        result=[
            JobResult(
                protein_id=item.protein_id,
                soluble_probability=item.soluble_probability
            )
            for item in job.results
        ],
        experiment_id=job.experiment.id
    )

class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            job_id = JobId(job_id)
            job: SolubilityJob = SolubilityJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            return map_job_to_response(job)
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e

            raise NoLabsException(ErrorCodes.unknown_exception) from e


class RunJobFeature:
    _api: DefaultApi

    def __init__(self, api: DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: SolubilityJob = SolubilityJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        result: List[Tuple[Protein, float]] = []

        for protein in job.proteins:
            response = self._api.run_solubility_run_post(
                run_solubility_prediction_request=RunSolubilityPredictionRequest(
                    amino_acid_sequence=protein.get_amino_acid_sequence(),
                    job_id=str(job_id.value)
                )
            )

            if response.errors:
                raise NoLabsException(ErrorCodes.amino_acid_localisation_run_error)

            soluble_probability = response.soluble_probability

            result.append((protein, soluble_probability))

        job.set_result(result)
        job.save(cascade=True)

        for protein, prob in result:
            protein.set_solubility_probability(SolubleProbability(prob))
            protein.save(cascade=True)

        return map_job_to_response(job)


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """
    async def handle(self, request: SetupJobRequest) -> JobResponse:
        assert request

        job_id = JobId(request.job_id if request.job_id else generate_uuid())
        job_name = JobName(request.job_name if request.job_name else 'New solubility job')

        jobs: SolubilityJob = SolubilityJob.objects(Q(id=job_id.value) | Q(name=job_name.value))

        if not jobs:
            if not request.experiment_id:
                raise NoLabsException(ErrorCodes.invalid_experiment_id)

            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job = SolubilityJob(
                id=job_id,
                name=job_name,
                experiment=experiment
            )
        else:
            job = jobs[0]

        proteins: List[Protein] = []
        for protein_id in request.protein_ids:
            protein = Protein.objects.with_id(protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            proteins.append(protein)

        job.set_input(proteins=proteins)
        job.save(cascade=True)

        return map_job_to_response(job)


class GetJobStatusFeature:
    def __init__(self, api: DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        job_id = JobId(job_id)

        job: SolubilityJob = SolubilityJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        response = self._api.is_job_running_job_job_id_is_running_get(job_id=str(job_id))

        return GetJobStatusResponse(
            running=response.is_running,
            result_valid=job.result_valid()
        )