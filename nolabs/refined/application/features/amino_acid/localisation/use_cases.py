__all__ = [
    'GetJobFeature',
    'RunLocalisationFeature'
]

import uuid
from typing import List
from uuid import UUID

from localisation_microservice import Configuration, DefaultApi, ApiClient, RunLocalisationPredictionRequest

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.controllers.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.controllers.amino_acid.localisation.api_models import GetJobResponse, \
    JobPropertiesResponse, \
    AminoAcidResponse, JobFastaPropertyResponse, RunJobResponse, UpdateJobRequest, GetJobMetadataResponse
from nolabs.refined.application.features.amino_acid.services import get_input_amino_acids
from nolabs.refined.domain.models import LocalisationJob
from nolabs.refined.domain.models.common import JobId, AminoAcid, Experiment, LocalisationProbability, AminoAcidId, \
    AminoAcidName, ExperimentId, \
    JobName
from nolabs.refined.infrastructure.settings import Settings


def map_to_amino_acid_response(amino_acid: AminoAcid) -> AminoAcidResponse:
    return AminoAcidResponse(
        sequence=amino_acid.get_content(),
        name=str(amino_acid.name),
        cytosolic=amino_acid.localisation.cytosolic,
        mitochondrial=amino_acid.localisation.mitochondrial,
        extracellular=amino_acid.localisation.extracellular,
        other=amino_acid.localisation.other,
        nuclear=amino_acid.localisation.nuclear
    )


# region Jobs

class UpdateJobFeature:
    async def handle(self, job_id: UUID, request: UpdateJobRequest):
        job_id = JobId(job_id)
        job: LocalisationJob = LocalisationJob.objects.get(job_id=job_id)

        updated = False
        if job.name != request.job_name:
            updated = True
            job.name = JobName(request.job_name)

        if updated:
            job.save()


class DeleteJobFeature:
    async def handle(self, job_id: UUID):
        job_id = JobId(job_id)
        job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)
        if not job:
            return
        job.delete()


class GetJobsMetadataFeature:
    async def handle(self, experiment_id: UUID) -> List[GetJobMetadataResponse]:
        assert experiment_id

        experiment_id = ExperimentId(experiment_id)
        experiment = Experiment.objects.with_id(experiment_id.value)

        if not experiment:
            raise NoLabsException('Experiment not found', ErrorCodes.experiment_not_found)

        result: List[GetJobMetadataResponse] = []
        job: LocalisationJob
        for job in LocalisationJob.objects.filter(experiment=experiment):
            result.append(
                GetJobMetadataResponse(
                    job_id=job.id.value,
                    job_name=str(job.name)
                )
            )

        return result


class GetJobFeature:
    async def handle(self, job_id: UUID) -> GetJobResponse:
        job_id = JobId(job_id)
        job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)

        if not job:
            NoLabsException.throw(ErrorCodes.job_not_found)

        amino_acids = job.amino_acids

        return GetJobResponse(
            job_id=job_id.value,
            job_name=str(job.name),
            amino_acids=[
                map_to_amino_acid_response(amino_acid)
                for amino_acid in job.amino_acids
            ],
            properties=JobPropertiesResponse(
                fastas=[
                    JobFastaPropertyResponse(
                        filename=f.name.fasta_name,
                        content=f.get_content()
                    ) for f in amino_acids
                ]
            )
        )


# endregion


class RunLocalisationFeature:
    _settings: Settings

    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self, request: RunAminoAcidRequest) -> RunJobResponse:
        assert request
        assert request.experiment_id

        if not request.fastas:
            raise NoLabsException('You did not specify amino acids', ErrorCodes.invalid_job_input)

        experiment_id = ExperimentId(request.experiment_id)
        experiment: Experiment = Experiment.objects.with_id(experiment_id.value)

        if not experiment:
            raise NoLabsException('Experiment not found', ErrorCodes.experiment_not_found)

        if not request.job_id:
            job = LocalisationJob(
                id=JobId(uuid.uuid4()),
                name=JobName('New job'),
                experiment=experiment
            )
            job.save()
            job_id = job.id
        else:
            job_id = JobId(request.job_id)
            job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)

        try:
            input_amino_acids = await get_input_amino_acids(experiment=experiment, request=request)

            job.set_amino_acids(input_amino_acids)
            job.save(cascade=True)

            configuration = Configuration(
                host=self._settings.localisation.microservice
            )

            job.clear_result()

            with ApiClient(configuration=configuration) as client:
                api_instance = DefaultApi(client)

                for amino_acid in job.amino_acids:
                    sequence = amino_acid.get_content()
                    result = api_instance.run_localisation_prediction_run_localisation_prediction_post(
                        run_localisation_prediction_request=RunLocalisationPredictionRequest(
                            amino_acid_sequence=sequence
                        )
                    )

                    if result.errors:
                        raise NoLabsException(result.errors, ErrorCodes.amino_acid_localisation_run_error)

                    localisation_probability = LocalisationProbability(
                        cytosolic=result.cytosolic_proteins,
                        mitochondrial=result.mitochondrial_proteins,
                        nuclear=result.nuclear_proteins,
                        other=result.other_proteins,
                        extracellular=result.extracellular_secreted_proteins
                    )

                    job.set_result(amino_acid, localisation=localisation_probability)
                    amino_acid.set_localisation_probability(localisation=localisation_probability)
                    amino_acid.save()

            job.save(cascade=True)

            results: List[AminoAcidResponse] = []

            for amino_acid in job.amino_acids:
                results.append(map_to_amino_acid_response(amino_acid))

            return RunJobResponse(job_id=job_id.value,
                                  amino_acids=results)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException('Unknown localisation error', ErrorCodes.unknown_localisation_error)
            raise e

