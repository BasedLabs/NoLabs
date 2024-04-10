from typing import Optional

from reinvent_microservice import Configuration, ApiClient, ReinventApi

from nolabs.api_models.small_molecules_design import GetExperimentResponse, ExperimentPropertiesResponse, \
    GetExperimentStatusResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.infrastructure.settings import Settings
from nolabs.modules.small_molecules_design.services.file_management import FileManagement


class GetExperimentFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._settings = settings
        self._fm = file_management

    async def handle(self, experiment_id: str) -> Optional[GetExperimentResponse]:
        if not self._fm.experiment_exists(experiment_id):
            return None

        configuration = Configuration(
            host=self._settings.reinvent_host,
        )

        experiment = self._fm.get_metadata(ExperimentId(experiment_id))

        params = self._fm.get_params(ExperimentId(experiment_id))
        target = self._fm.get_pdb(ExperimentId(experiment_id))

        if not target:
            properties = ExperimentPropertiesResponse(
                pdb_file=None,
                pdb_file_name=None,
                center_x=0.0,
                center_y=0.0,
                center_z=0.0,
                size_x=5.0,
                size_y=5.0,
                size_z=5.0,
                batch_size=128,
                minscore=0.4,
                epochs=50
            )
        else:
            properties = ExperimentPropertiesResponse(
                pdb_file=target.read_string(),
                pdb_file_name=target.name,
                center_x=params.center_x,
                center_y=params.center_y,
                center_z=params.center_z,
                size_x=params.size_x,
                size_y=params.size_y,
                size_z=params.size_z,
                batch_size=params.batch_size,
                minscore=params.minscore,
                epochs=params.epochs,
            )

        with ApiClient(configuration=configuration) as client:
            api_instance = ReinventApi(client)

            config_result = api_instance.get_config_api_reinvent_reinvent_config_id_get(experiment_id)

            if not config_result:
                return GetExperimentResponse(
                    experiment_id=experiment_id,
                    experiment_name=experiment.name.value,
                    created_at=experiment.date,
                    status=GetExperimentStatusResponse(running=False, sampling_allowed=False),
                    properties=properties
                )

            config = config_result.actual_instance

            running = config.running if config_result else False
            sampling_allowed = config.sampling_allowed if config_result else False

            return GetExperimentResponse(
                experiment_id=experiment_id,
                experiment_name=experiment.name.value,
                created_at=experiment.date,
                status=GetExperimentStatusResponse(running=running, sampling_allowed=sampling_allowed),
                properties=properties
            )
