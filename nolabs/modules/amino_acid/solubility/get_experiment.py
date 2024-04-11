from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.amino_acid.solubility import GetExperimentResponse, ExperimentPropertiesResponse, \
    ExperimentFastaPropertyResponse
from nolabs.modules.amino_acid.solubility.services.file_management import FileManagement


class GetExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    async def handle(self, id: str) -> GetExperimentResponse:
        assert id

        experiment_id = ExperimentId(id)

        if not self._file_management.metadata_exists(experiment_id) or not self._file_management.properties_exists(
                experiment_id):
            raise NoLabsException(['This experiment not found'], ErrorCodes.experiment_not_found)

        data = self._file_management.get_result(experiment_id)
        properties = await self._file_management.get_properties(experiment_id)
        metadata = self._file_management.get_metadata(experiment_id)
        return GetExperimentResponse(
            experiment_id=experiment_id.value,
            experiment_name=metadata.name.value,
            amino_acids=data,
            properties=ExperimentPropertiesResponse(
                amino_acid_sequence=properties.amino_acid_sequence,
                fastas=[
                    ExperimentFastaPropertyResponse(
                        filename=f.filename,
                        content=(await f.read()).decode('utf-8')
                    ) for f in properties.fastas
                ]
            )
        )
