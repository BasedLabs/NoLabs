from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.conformations import GetExperimentResponse, \
    ExperimentPropertiesResponse, IntegratorsRequest
from nolabs.modules.conformations.services.file_management import FileManagement
from nolabs.exceptions import NoLabsException, ErrorCodes


class GetExperimentFeature:
    def __init__(self, file_management: FileManagement):
        self._file_management = file_management

    async def handle(self, id: str) -> GetExperimentResponse:
        assert id

        experiment_id = ExperimentId(id)

        if not self._file_management.metadata_exists(experiment_id) or not self._file_management.properties_exists(
                experiment_id):
            raise NoLabsException(messages=["Experiment does not exist"], error_code=ErrorCodes.experiment_not_found)

        metadata = self._file_management.get_metadata(experiment_id)
        properties = self._file_management.get_properties(experiment_id)
        data, timeline = self._file_management.get_result(experiment_id)
        return GetExperimentResponse(
            experiment_id=experiment_id.value,
            experiment_name=metadata.name.value,
            pdb_file=data,
            timeline=timeline,
            properties=ExperimentPropertiesResponse(
                pdb_file=(await properties.pdb_file.read()).decode('utf-8'),
                pdb_file_name=properties.pdb_file.filename,
                total_frames=properties.total_frames,
                temperature_k=properties.temperature_k,
                take_frame_every=properties.take_frame_every,
                step_size=properties.step_size,
                replace_non_standard_residues=properties.replace_non_standard_residues,
                add_missing_atoms=properties.add_missing_atoms,
                add_missing_hydrogens=properties.add_missing_hydrogens,
                friction_coeff=properties.friction_coeff,
                ignore_missing_atoms=properties.ignore_missing_atoms,
                integrator=IntegratorsRequest(properties.integrator)
            )
        )
