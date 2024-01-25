from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.protein_design import GetExperimentResponse, \
    ExperimentPropertiesResponse
from nolabs.features.protein_design.services.file_management import FileManagement
from nolabs.exceptions import NoLabsException, ErrorCodes


class GetExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> GetExperimentResponse:
        assert id

        experiment_id = ExperimentId(id)

        if not self._file_management.experiment_exists(experiment_id):
            raise NoLabsException(messages=["Experiment does not exist"], error_code=ErrorCodes.experiment_id_not_found)

        metadata = self._file_management.get_experiment_metadata(experiment_id)
        properties = self._file_management.experiment_properties(experiment_id)
        data = self._file_management.get_experiment_data(experiment_id)
        return GetExperimentResponse(
            experiment_id=experiment_id.value,
            experiment_name=metadata.name.value,
            pdb_files=data,
            properties=ExperimentPropertiesResponse(
                pdb_file=properties.input_file_content,
                contig=properties.contig,
                number_of_designs=properties.number_of_designs,
                hotspots=properties.hotspots,
                timesteps=properties.timesteps,
                pdb_file_name=properties.input_file_name
            )
        )
