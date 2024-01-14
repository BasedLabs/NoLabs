from typing import List, Tuple

from solubility_microservice import (RunSolubilityPredictionRequest,
                                     ApiClient, DefaultApi, Configuration)

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.amino_acid import AminoAcid
from nolabs.domain.solubility import SolubilityProbability
from nolabs.infrastructure.settings import Settings
from nolabs.api_models.solubility import RunSolubilityResponse, RunSolubilityRequest, AminoAcidResponse
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.features.solubility.services.file_management import FileManagement
from nolabs.utils import generate_uuid
from nolabs.utils.fasta import FastaReader


class RunSolubilityFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._settings = settings
        self._file_management = file_management

    async def handle(self, request: RunSolubilityRequest) -> RunSolubilityResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id) if request.experiment_id else ExperimentId(generate_uuid())

        configuration = Configuration(
            host=self._settings.solubility_host,
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
                    for amino_acid in fasta_reader.get_ids_and_sequences(fasta_content.decode('utf-8')):
                        if amino_acid.name not in [aa.name for aa in amino_acids]:
                            amino_acids.append(amino_acid)

            results: List[Tuple[AminoAcid, SolubilityProbability]] = []
            if not amino_acids:
                raise NoLabsException('No amino acids', ErrorCodes.no_amino_acids)
            for amino_acid in amino_acids:
                result = api_instance.run_solubility_run_solubility_prediction_post(
                    run_solubility_prediction_request=RunSolubilityPredictionRequest(
                        amino_acid_sequence=amino_acid.sequence
                    )
                )
                if result.errors:
                    raise NoLabsException(','.join(result.errors), ErrorCodes.amino_acid_solubility_run_error)
                results.append((amino_acid, SolubilityProbability(result.soluble_probability)))  # type: ignore
            self._file_management.update_metadata(experiment_id=experiment_id,
                                                  experiment_name=ExperimentName(value=request.experiment_name),
                                                  amino_acids=amino_acids)
            self._file_management.save_experiment(
                experiment_id=experiment_id,
                data=results
            )
            return RunSolubilityResponse(experiment_id=experiment_id.value,
                                         amino_acids=[AminoAcidResponse(
                                             sequence=amino_acid.sequence,
                                             name=amino_acid.name,
                                             soluble_probability=prob.value
                                         ) for amino_acid, prob in results])



