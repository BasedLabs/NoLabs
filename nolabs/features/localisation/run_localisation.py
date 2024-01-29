from typing import List, Tuple

from localisation_microservice import (RunLocalisationPredictionRequest,
                                       ApiClient, DefaultApi, Configuration)

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.amino_acid import AminoAcid
from nolabs.infrastructure.settings import Settings
from nolabs.api_models.localisation import RunLocalisationResponse, RunLocalisationRequest, AminoAcidResponse
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.features.localisation.services.file_management import FileManagement
from nolabs.utils.fasta import FastaReader


class RunLocalisationFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._settings = settings
        self._file_management = file_management

    async def handle(self, request: RunLocalisationRequest) -> RunLocalisationResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        self._file_management.ensure_experiment_folder_exists(experiment_id=experiment_id)
        self._file_management.cleanup_experiment(experiment_id=experiment_id)
        await self._file_management.set_metadata(experiment_id=experiment_id, experiment_name=ExperimentName(request.experiment_name))
        await self._file_management.set_properties(experiment_id=experiment_id, request=request)

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

            self._file_management.set_metadata(experiment_id=experiment_id,
                                               experiment_name=ExperimentName(value=request.experiment_name))
            self._file_management.set_result(experiment_id=experiment_id, data=results)
            return RunLocalisationResponse(experiment_id=experiment_id.value,
                                           experiment_name=request.experiment_name,
                                           amino_acids=results)
