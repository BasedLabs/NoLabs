from reinvent_microservice import Configuration, ApiClient, ReinventApi

from nolabs.infrastructure.settings import Settings


class StartSamplingExperimentFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self, experiment_id: str):
        configuration = Configuration(
            host=self._settings.reinvent_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = ReinventApi(client)
            api_instance.sampling_api_reinvent_config_id_start_sampling_post(config_id=experiment_id)
