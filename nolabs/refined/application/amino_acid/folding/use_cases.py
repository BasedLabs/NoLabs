__all__ = [
    'GetFoldingJobFeature',
    'RunFoldingFeature'
]

import uuid
from typing import List
from uuid import UUID

from esmfold_microservice import DefaultApi, RunEsmFoldPredictionRequest

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.amino_acid.folding.api_models import GetJobResponse, \
    JobPropertiesResponse, \
    AminoAcidResponse, JobFastaPropertyResponse, RunJobResponse
from nolabs.refined.application.amino_acid.services import get_input_proteins
from nolabs.refined.domain.models import LocalisationJob
from nolabs.refined.domain.models.common import JobId, Experiment, ExperimentId, \
    JobName, Protein
from nolabs.refined.domain.models.folding import FoldingJob, FoldingBackendEnum
from nolabs.refined.infrastructure.settings import Settings


def map_to_amino_acid_response(protein: Protein) -> AminoAcidResponse:
    return AminoAcidResponse(
        sequence=protein.get_fasta(),
        name=str(protein.name),
        pdb_file_name=protein.name.pdb_name,
        pdb_file=protein.get_pdb()
    )


class GetFoldingJobFeature:
    async def handle(self, job_id: UUID) -> GetJobResponse:
        job_id = JobId(job_id)
        job: LocalisationJob = LocalisationJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        proteins = job.proteins

        return GetJobResponse(
            job_id=job_id.value,
            job_name=str(job.name),
            amino_acids=[
                map_to_amino_acid_response(protein)
                for protein in job.proteins
            ],
            properties=JobPropertiesResponse(
                fastas=[
                    JobFastaPropertyResponse(
                        filename=f.name.fasta_name,
                        content=f.get_fasta()
                    ) for f in proteins
                ]
            )
        )


class RunFoldingFeature:
    _settings: Settings
    _api: DefaultApi

    def __init__(self, api: DefaultApi, settings: Settings):
        self._settings = settings
        self._api = api

    async def handle(self, request: RunAminoAcidRequest) -> RunJobResponse:
        assert request
        assert request.experiment_id

        if not request.fastas:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        experiment_id = ExperimentId(request.experiment_id)
        experiment: Experiment = Experiment.objects.with_id(experiment_id.value)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        if not request.job_id:
            job = FoldingJob(
                id=JobId(uuid.uuid4()),
                name=JobName('New job'),
                experiment=experiment
            )
            job.save()
            job_id = JobId(job.id)
        else:
            job_id = JobId(request.job_id)
            job: FoldingJob = FoldingJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        try:
            input_proteins = await get_input_proteins(experiment=experiment, request=request)

            job.set_proteins(input_proteins)
            job.save(cascade=True)

            job.clear_result()

            for protein in job.proteins:
                sequence = protein.get_fasta()
                result = self._api.predict_run_folding_post(
                    run_esm_fold_prediction_request=RunEsmFoldPredictionRequest(
                        protein_sequence=sequence
                    ), _request_timeout=(1000.0, 1000.0))

                if result.errors:
                    raise NoLabsException(ErrorCodes.amino_acid_localisation_run_error)

                job.set_result(protein, FoldingBackendEnum.esmfold, result.pdb_content)
                protein.set_pdb(pdb_content=result.pdb_content)
                protein.save()

            job.save(cascade=True)

            results: List[AminoAcidResponse] = []

            for protein in job.proteins:
                results.append(map_to_amino_acid_response(protein))

            return RunJobResponse(job_id=job_id.value,
                                  amino_acids=results)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_folding_error) from e
            raise e
