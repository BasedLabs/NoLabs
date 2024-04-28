__all__ = [
    'GetJobFeature',
    'RunLocalisationFeature'
]

from typing import List
from uuid import UUID

from localisation_microservice import Configuration, DefaultApi, ApiClient, RunLocalisationPredictionRequest

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.controllers.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.controllers.amino_acid.localisation.api_models import GetJobResponse, \
    JobPropertiesResponse, \
    AminoAcidResponse, JobFastaPropertyResponse, RunJobResponse, UpdateJobRequest, GetJobMetadataResponse
from nolabs.refined.domain.models import LocalisationJob, LocalisationProbability
from nolabs.refined.domain.models.common import JobId, AminoAcid, JobType, AminoAcidId, AminoAcidName, ExperimentId, \
    JobName
from nolabs.refined.domain.models.experiment import Experiment
from nolabs.refined.domain.repository import Repository
from nolabs.refined.infrastructure.settings import Settings
from nolabs.utils.fasta import FastaReader


def map_probability_to_amino_acid_response(probability: LocalisationProbability) -> AminoAcidResponse:
    return AminoAcidResponse(
        sequence=probability.amino_acid.content.read(),
        name=probability.amino_acid.name.value,
        cytosolic_proteins=probability.cytosolic_proteins,
        mitochondrial_proteins=probability.mitochondrial_proteins,
        extracellular_secreted_proteins=probability.extracellular_secreted_proteins,
        other_proteins=probability.other_proteins,
        nuclear_proteins=probability.nuclear_proteins
    )

# region Jobs

class UpdateJobFeature:
    _repository: Repository

    def __init__(self, repository: Repository):
        self._repository = repository

    async def handle(self, job_id: UUID, request: UpdateJobRequest):
        job_id = JobId(job_id)
        job: LocalisationJob = self._repository.localisation_jobs.get(job_id=job_id)

        updated = False
        if job.name != request.job_name:
            updated = True
            job.name = JobName(request.job_name)

        if updated:
            job.save()


class DeleteJobFeature:
    _repository: Repository

    def __init__(self, repository: Repository):
        self._repository = repository

    async def handle(self, job_id: UUID):
        job_id = JobId(job_id)
        self._repository.localisation_jobs.delete(job_id)


class GetJobsMetadataFeature:
    _repository: Repository

    def __init__(self, repository: Repository):
        self._repository = repository

    async def handle(self, experiment_id: UUID) -> List[GetJobMetadataResponse]:
        experiment_id = ExperimentId(experiment_id)
        experiment = self._repository.experiments(experiment_id)

        if not experiment:
            NoLabsException.throw(ErrorCodes.experiment_not_found)

        result: List[GetJobMetadataResponse] = []
        job: LocalisationJob
        for job in self._repository.localisation_jobs.filter(experiment=experiment):
            result.append(
                GetJobMetadataResponse(
                    job_id=job.id.value,
                    job_name=str(job.name)
                )
            )

        return result


class GetJobFeature:
    _repository: Repository

    def __init__(self, repository: Repository):
        self._repository = repository

    async def handle(self, job_id: UUID) -> GetJobResponse:
        job_id = JobId(job_id)
        job = self._repository.localisation_jobs.get(id=job_id)

        if not job:
            NoLabsException.throw(ErrorCodes.job_not_found)

        amino_acids = job.amino_acids

        return GetJobResponse(
            job_id=job_id.value,
            job_name=str(job.name),
            amino_acids=[
                map_probability_to_amino_acid_response(probability)
                for probability in job.result
            ],
            properties=JobPropertiesResponse(
                fastas=[
                    JobFastaPropertyResponse(
                        filename=f.name.fasta_name,
                        content=f.content.read()
                    ) for f in amino_acids
                ]
            )
        )

# endregion


class RunLocalisationFeature:
    _repository: Repository
    _settings: Settings

    def __init__(self, repository: Repository, settings: Settings):
        self._settings = settings
        self._repository = repository

    async def handle(self, request: RunAminoAcidRequest) -> RunJobResponse:
        assert request

        job_id = JobId(request.job_id)
        experiment_id = ExperimentId(request.experiment_id)
        experiment: Experiment = self._repository.experiments.get(id=experiment_id)

        if not experiment:
            raise NoLabsException('Experiment not found', ErrorCodes.experiment_not_found)

        job: LocalisationJob = self._repository.localisation_jobs.get(id=job_id)
        if not job:
            job = LocalisationJob(
                id=job_id,
                name=JobName('New job'),
                experiment=experiment
            )
            self._repository.localisation_jobs.insert(job)

        try:
            input_amino_acids = await self._get_input_amino_acids(experiment=experiment, request=request)

            job.set_amino_acids(input_amino_acids)
            experiment.save()

            output = await self._predictor(input_amino_acids)

            job.set_result(output)
            experiment.save()

            results: List[AminoAcidResponse] = []

            for probability in job.result:
                results.append(map_probability_to_amino_acid_response(probability))

            return RunJobResponse(job_id=job_id.value,
                                  amino_acids=results)
        except Exception as e:
            if not isinstance(e, NoLabsException):
                raise NoLabsException('Unknown localisation error', ErrorCodes.unknown_localisation_error)
            raise e

    async def _predictor(self, inputs: List[AminoAcid]) -> List[LocalisationProbability]:
        configuration = Configuration(
            host=self._settings.localisation.microservice
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)

            output: List[LocalisationProbability] = []
            for amino_acid in inputs:
                sequence = amino_acid.content.value.decode('utf-8')
                result = api_instance.run_localisation_prediction_run_localisation_prediction_post(
                    run_localisation_prediction_request=RunLocalisationPredictionRequest(
                        amino_acid_sequence=sequence
                    )
                )

                if result.errors:
                    raise NoLabsException(result.errors, ErrorCodes.amino_acid_localisation_run_error)

                output.append(
                    LocalisationProbability(
                        amino_acid=amino_acid,
                        cytosolic_proteins=result.cytosolic_proteins,
                        mitochondrial_proteins=result.mitochondrial_proteins,
                        nuclear_proteins=result.nuclear_proteins,
                        other_proteins=result.other_proteins,
                        extracellular_secreted_proteins=result.extracellular_secreted_proteins
                    )
                )

        return output

    async def _get_input_amino_acids(self, experiment: Experiment, request: RunAminoAcidRequest) -> List[AminoAcid]:
        input_amino_acids: List[AminoAcid] = []
        fasta_reader = FastaReader()
        for fasta in request.fastas:
            fasta_content = await fasta.read()
            for amino_acid in fasta_reader.get_ids2seqs(fasta_content.decode('utf-8')):
                existing_amino_acid = self._repository.amino_acids.get(name=AminoAcidName(amino_acid.name))
                if existing_amino_acid:
                    existing_amino_acid.content.put(amino_acid.sequence)
                else:
                    existing_amino_acid = AminoAcid(
                        id=AminoAcidId(UUID()),
                        name=AminoAcidName(amino_acid.name),
                        content=amino_acid.sequence,
                        experiment=experiment
                    )
                existing_amino_acid.save()
                input_amino_acids.append(existing_amino_acid)
        return input_amino_acids
