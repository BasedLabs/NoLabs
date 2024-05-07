__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature'
]

from typing import List, Tuple
from uuid import UUID

from solubility_microservice import DefaultApi, RunSolubilityPredictionRequest

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.solubility.api_models import JobResponse, JobResult, SetupJobRequest
from nolabs.refined.domain.models.common import JobId, Experiment, JobName, Protein, SolubleProbability
from nolabs.refined.domain.models.solubility import SolubilityJob
from nolabs.utils import generate_uuid


def map_job_to_response(job: SolubilityJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        proteins=[p.iid.value for p in job.proteins],
        result=[
            JobResult(
                protein_id=item.protein_id,
                soluble_probability=item.soluble_probability
            )
            for item in job.results
        ]
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
        try:
            assert job_id

            job_id = JobId(job_id)
            job: SolubilityJob = SolubilityJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            result: List[Tuple[Protein, float]] = []

            for protein in job.proteins:
                response = self._api.run_solubility_run_solubility_prediction_post(
                    run_solubility_prediction_request=RunSolubilityPredictionRequest(
                        amino_acid_sequence=protein.get_fasta()
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
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_exception) from e
            raise e


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """
    async def handle(self, request: SetupJobRequest) -> JobResponse:
        try:
            assert request

            job_id = JobId(request.job_id if request.job_id else generate_uuid())
            job_name = JobName(request.job_name if request.job_name else 'New solubility job')

            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job: SolubilityJob = SolubilityJob.objects.with_id(job_id.value)

            if not job:
                job = SolubilityJob(
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

            job.set_proteins(proteins=proteins)
            job.save(cascade=True)

            return map_job_to_response(job)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_solubility_error) from e
            raise e
