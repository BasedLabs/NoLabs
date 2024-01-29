from typing import List

from esmfold_microservice import (ApiClient, DefaultApi, Configuration, RunEsmFoldPredictionResponse,
                                  RunEsmFoldPredictionRequest)
from slugify import slugify

from nolabs.api_models.folding import RunFoldingRequest, RunFoldingResponse, AminoAcidResponse
from nolabs.domain.amino_acid import AminoAcid
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.features.folding.services.file_management import FileManagement
from nolabs.infrastructure.settings import Settings
from nolabs.utils.fasta import FastaReader


class RunFoldingFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._settings = settings
        self._file_management = file_management

    async def handle(self, request: RunFoldingRequest) -> RunFoldingResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        self._file_management.ensure_experiment_folder_exists(experiment_id=experiment_id)
        self._file_management.cleanup_experiment(experiment_id=experiment_id)
        await self._file_management.set_metadata(experiment_id=experiment_id,
                                                 experiment_name=ExperimentName(request.experiment_name))
        await self._file_management.set_properties(experiment_id=experiment_id, request=request)

        configuration = Configuration(
            host=self._settings.folding_host
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
                result = api_instance.predict_run_folding_post(
                    run_esm_fold_prediction_request=RunEsmFoldPredictionRequest(
                        protein_sequence=amino_acid.sequence
                    )
                )
                results.append(AminoAcidResponse(
                    sequence=amino_acid.sequence,
                    name=amino_acid.name,
                    pdb_file=result.pdb_content,
                    pdb_file_name=f'{slugify(amino_acid.name)}.pdb'
                ))  # type: ignore

            self._file_management.set_result(experiment_id=experiment_id, data=results)
            return RunFoldingResponse(experiment_id=experiment_id.value,
                                           experiment_name=request.experiment_name,
                                           amino_acids=results)
