from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.protein_design import GetExperimentResponse, ExperimentMetadataResponse
from nolabs.features.protein_design.services.file_management import FileManagement
from nolabs.exceptions import NoLabsException, ErrorCodes


class GetExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    def handle(self, id: str) -> GetExperimentResponse:
        assert id

        experiment_id = ExperimentId(id)

        if not self._file_management.experiment_exists(experiment_id):
            raise NoLabsException(message="Experiment does not exist", error_code=ErrorCodes.experiment_id_not_found)

        metadata = self._file_management.get_experiment_metadata(experiment_id)
        input_file = self._file_management.get_input_pdb_file(experiment_id)
        data = self._file_management.get_experiment_data(experiment_id)
        return GetExperimentResponse(
            experiment_id=experiment_id.value,
            experiment_name=metadata.name.value,
            pdbs_content=data,
            pdb_file=input_file,
            pdb_file_name=metadata.properties['pdb_file_name'],

        )
