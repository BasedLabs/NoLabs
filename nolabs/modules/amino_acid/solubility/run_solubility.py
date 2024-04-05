from typing import List

from solubility_microservice import (RunSolubilityPredictionRequest,
                                     ApiClient, DefaultApi, Configuration)

from nolabs.api_models.amino_acid.common_models import RunAminoAcidRequest
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.amino_acid import AminoAcid
from nolabs.modules.amino_acid.run_aa_inference_feature_base import RunAminoAcidInferenceFeature
from nolabs.infrastructure.settings import Settings
from nolabs.api_models.amino_acid.solubility import RunSolubilityResponse, AminoAcidResponse
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.modules.amino_acid.solubility.services.file_management import FileManagement
from nolabs.utils.fasta import FastaReader


class RunSolubilityFeature(RunAminoAcidInferenceFeature[FileManagement]):
    def __init__(self, settings: Settings, file_management: FileManagement):
        super().__init__(file_management)
        self._settings = settings

    async def handle(self, request: RunAminoAcidRequest) -> RunSolubilityResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        await self._setup_experiment(experiment_id, request)

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
                    for amino_acid in fasta_reader.get_ids2seqs(fasta_content.decode('utf-8')):
                        if amino_acid.name not in [aa.name for aa in amino_acids]:
                            amino_acids.append(amino_acid)

            results: List[AminoAcidResponse] = []
            if not amino_acids:
                raise NoLabsException(['No amino acids'], ErrorCodes.no_amino_acids)
            for amino_acid in amino_acids:
                result = api_instance.run_solubility_run_solubility_prediction_post(
                    run_solubility_prediction_request=RunSolubilityPredictionRequest(
                        amino_acid_sequence=amino_acid.sequence
                    )
                )
                if result.errors:
                    raise NoLabsException(result.errors, ErrorCodes.amino_acid_solubility_run_error)
                results.append(AminoAcidResponse(
                    sequence=amino_acid.sequence,
                    name=amino_acid.name,
                    soluble_probability=result.soluble_probability
                ))  # type: ignore
            self._file_management.set_result(
                experiment_id=experiment_id,
                data=results
            )
            metadata = self._file_management.get_metadata(experiment_id)
            return RunSolubilityResponse(experiment_id=experiment_id.value,
                                         experiment_name=metadata.name.value,
                                         amino_acids=results, errors=[])
