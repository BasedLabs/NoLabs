from reinvent_microservice import Configuration, ApiClient, ReinventApi

from nolabs.api_models.small_molecules_design import GetExperimentStatusResponse
from nolabs.infrastructure.settings import Settings
from nolabs.modules.small_molecules_design.services.file_management import FileManagement


class GetExperimentStatusFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._settings = settings
        self._fm = file_management

    async def handle(self, experiment_id: str) -> GetExperimentStatusResponse:
        if not self._fm.experiment_exists(experiment_id):
            return GetExperimentStatusResponse(running=False, sampling_allowed=False)

        configuration = Configuration(
            host=self._settings.reinvent_host,
        )

        with ApiClient(configuration=configuration) as client:
            api_instance = ReinventApi(client)

            config_result = api_instance.get_config_api_reinvent_reinvent_config_id_get(experiment_id)

            if not config_result:
                return GetExperimentStatusResponse(running=False, sampling_allowed=False)

            config = config_result.actual_instance

            return GetExperimentStatusResponse(running=config.running, sampling_allowed=config.sampling_allowed)
