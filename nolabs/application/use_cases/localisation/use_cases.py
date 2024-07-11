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
from nolabs.application.use_cases.localisation.api_models import JobResponse, JobResult, SetupJobRequest, \
    GetJobStatusResponse
from nolabs.domain.models.localisation import LocalisationJob
from nolabs.domain.models.common import JobId, Experiment, LocalisationProbability, JobName, Protein
from nolabs.utils import generate_uuid


def map_job_to_response(job: LocalisationJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        protein_ids=[p.iid.value for p in job.proteins],
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
        ],
        experiment_id=job.experiment.id
    )


class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        job_id = JobId(job_id)
        job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return map_job_to_response(job)


class RunJobFeature:
    _api: DefaultApi

    def __init__(self, api: DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        result: List[Tuple[Protein, LocalisationProbability]] = []

        try:
            job.started()
            await job.save()

            for protein in job.proteins:
                response = self._api.predict_run_post(
                    run_localisation_prediction_request=RunLocalisationPredictionRequest(
                        amino_acid_sequence=protein.get_amino_acid_sequence(),
                        job_id=str(job_id.value)
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
            await job.save(cascade=True)

            for protein, prob in result:
                protein.set_localisation_probability(prob)
                protein.save(cascade=True)
        finally:
            job.finished()
            await job.save()

        return map_job_to_response(job)


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """

    async def handle(self, request: SetupJobRequest) -> JobResponse:
        assert request

        job_id = JobId(request.job_id if request.job_id else generate_uuid())
        job_name = JobName(request.job_name if request.job_name else 'New localisation job')

        jobs: LocalisationJob = LocalisationJob.objects(Q(id=job_id.value) | Q(name=job_name.value))

        if not jobs:
            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job = LocalisationJob(
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
        await job.save(cascade=True)

        return map_job_to_response(job)


class GetJobStatusFeature:
    def __init__(self, api: DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> GetJobStatusResponse:
        job_id = JobId(job_id)

        job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        response = self._api.is_job_running_job_job_id_is_running_get(job_id=str(job_id))

        return GetJobStatusResponse(
            running=response.is_running,
            result_valid=job.result_valid()
        )
