from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.solubility import GetExperimentResponse, ExperimentMetadataResponse, AminoAcidResponse
from nolabs.features.solubility.services.file_management import FileManagement


class GetExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> GetExperimentResponse:
        assert id

        experiment_id = ExperimentId(id)

        if not self._file_management.experiment_exists(experiment_id):
            raise NoLabsException('This experiment not found in solubility', ErrorCodes.experiment_id_not_found)

        metadata = self._file_management.get_metadata(experiment_id)
        data = self._file_management.get_experiment_data(experiment_id)
        return GetExperimentResponse(
            metadata=ExperimentMetadataResponse(
                id=metadata.id.value,
                name=metadata.name.value,
                date=metadata.date
            ),
            amino_acids=[AminoAcidResponse(
                sequence=amino_acid.sequence,
                name=amino_acid.name,
                soluble_probability=prob.value
            ) for amino_acid, prob in data]
        )
