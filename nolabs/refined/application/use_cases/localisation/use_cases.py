__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature'
]

from typing import List, Tuple
from uuid import UUID

from localisation_microservice import DefaultApi, RunLocalisationPredictionRequest
from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.localisation.api_models import JobResponse, JobResult, SetupJobRequest
from nolabs.refined.domain.models.localisation import LocalisationJob
from nolabs.refined.domain.models.common import JobId, Experiment, LocalisationProbability, JobName, Protein
from nolabs.utils import generate_uuid


def map_job_to_response(job: LocalisationJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        proteins=[p.iid.value for p in job.proteins],
        result=[
            JobResult(
                protein_id=item.protein_id,
                cytosolic=item.cytosolic,
                mitochondrial=item.mitochondrial,
                nuclear=item.nuclear,
                other=item.other,
                extracellular=item.extracellular
            )
            for item in job.probabilities
        ]
    )

class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        try:
            job_id = JobId(job_id)
            job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)

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
            job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            result: List[Tuple[Protein, LocalisationProbability]] = []

            for protein in job.proteins:
                response = self._api.run_localisation_prediction_run_localisation_prediction_post(
                    run_localisation_prediction_request=RunLocalisationPredictionRequest(
                        amino_acid_sequence=protein.get_amino_acid_sequence()
                    )
                )

                if response.errors:
                    raise NoLabsException(ErrorCodes.amino_acid_localisation_run_error)

                localisation_probability = LocalisationProbability(
                    cytosolic=response.cytosolic_proteins,
                    mitochondrial=response.mitochondrial_proteins,
                    nuclear=response.nuclear_proteins,
                    other=response.other_proteins,
                    extracellular=response.extracellular_secreted_proteins
                )

                result.append((protein, localisation_probability))

            job.set_result(result)
            job.save(cascade=True)

            for protein, prob in result:
                protein.set_localisation_probability(prob)
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
            job_name = JobName(request.job_name if request.job_name else 'New localisation job')

            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job: LocalisationJob = LocalisationJob.objects(Q(id=job_id.value) | Q(name=job_name.value)).first()

            if not job:
                job = LocalisationJob(
                    id=job_id,
                    name=job_name,
                    experiment=experiment
                )

            proteins: List[Protein] = []
            for protein_id in request.proteins:
                protein = Protein.objects(id=protein_id, experiment=experiment).first()

                if not protein:
                    raise NoLabsException(ErrorCodes.protein_not_found)

                proteins.append(protein)

            job.set_input(proteins=proteins)
            job.save(cascade=True)

            return map_job_to_response(job)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_folding_error) from e
            raise e

