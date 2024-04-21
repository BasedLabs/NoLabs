__all__ = [
    'GetJobFeature',
    'RunLocalisationFeature'
]

from typing import List
from uuid import UUID

from localisation_microservice import Configuration, DefaultApi, ApiClient

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.infrastructure.settings import Settings
from nolabs.refined.application.features.amino_acid.localisation.api_models import GetJobResponse, JobPropertiesResponse, \
    AminoAcidResponse, JobFastaPropertyResponse

from nolabs.refined.application.features.amino_acid.common_api_models import RunAminoAcidRequest
from nolabs.refined.application.features.amino_acid.localisation.api_models import RunJobResponse
from nolabs.refined.domain.common.entities import JobId, AminoAcid
from nolabs.refined.domain.repository import Repository


class GetJobFeature:
    def __init__(self, repository: Repository):
        self._repository = repository

    async def handle(self, job_id: UUID) -> GetJobResponse:
        job_id = JobId(job_id)
        job = self._repository.get_localisation_job(job_id)

        if not job:
            NoLabsException.throw(ErrorCodes.job_not_found)

        amino_acids = self._repository.get_amino_acids_by_ids(job.amino_acid_ids)

        amino_acids_response = [
            AminoAcidResponse(
                sequence=str(amino_acids[probability.amino_acid_id].content),
                name=amino_acids[probability.amino_acid_id].name.value,
                cytosolic_proteins=probability.cytosolic_proteins,
                mitochondrial_proteins=probability.mitochondrial_proteins,
                extracellular_secreted_proteins=probability.extracellular_secreted_proteins,
                other_proteins=probability.other_proteins,
                nuclear_proteins=probability.nuclear_proteins
            )
            for probability in job.output
        ]

        return GetJobResponse(
            job_id=job_id.value,
            job_name=str(job.name),
            amino_acids=amino_acids_response,
            properties=JobPropertiesResponse(
                fastas=[
                    JobFastaPropertyResponse(
                        filename=f.name.fasta_name,
                        content=str(f.content)
                    ) for f in amino_acids.values()
                ]
            )
        )


class RunLocalisationFeature:
    _settings: Settings
    _repository: Repository

    def __init__(self, settings: Settings, repository: Repository):
        self._settings = settings
        self._repository = repository

    async def handle(self, request: RunAminoAcidRequest) -> RunJobResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        await self._setup_experiment(experiment_id=experiment_id, request=request)

        configuration = Configuration(
            host=self._settings.localisation_host
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            amino_acids: List[AminoAcid] = []

            if request.amino_acid_sequence:
                amino_acids.append(AminoAcid(name='unknown', sequence=request.amino_acid_sequence))

            if request.fastas:
                fasta_reader = FastaReader()
                for fasta in request.fastas:
                    fasta_content = await fasta.read()
                    for amino_acid in fasta_reader.get_ids2seqs(fasta_content.decode('utf-8')):
                        if amino_acid.name not in [aa.name for aa in amino_acids]:
                            amino_acids.append(amino_acid)

            results: List[AminoAcidResponse] = []

            if not amino_acids:
                raise NoLabsException(['No amino acids'], ErrorCodes.no_amino_acids)

            for amino_acid in amino_acids:
                result = api_instance.run_localisation_prediction_run_localisation_prediction_post(
                    run_localisation_prediction_request=RunLocalisationPredictionRequest(
                        amino_acid_sequence=amino_acid.sequence
                    )
                )
                if result.errors:
                    raise NoLabsException(result.errors, ErrorCodes.amino_acid_localisation_run_error)
                results.append(AminoAcidResponse(
                    sequence=amino_acid.sequence,
                    name=amino_acid.name,
                    cytosolic_proteins=result.cytosolic_proteins,
                    extracellular_secreted_proteins=result.extracellular_secreted_proteins,
                    nuclear_proteins=result.nuclear_proteins,
                    other_proteins=result.other_proteins,
                    mitochondrial_proteins=result.mitochondrial_proteins
                ))  # type: ignore

            self._file_management.set_result(experiment_id=experiment_id, data=results)
            return RunLocalisationResponse(experiment_id=experiment_id.value,
                                           experiment_name=request.experiment_name,
                                           amino_acids=results)
