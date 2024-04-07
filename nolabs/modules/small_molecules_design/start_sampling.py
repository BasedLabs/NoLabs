from reinvent_microservice import Configuration, ApiClient, DefaultApi

from nolabs.infrastructure.settings import Settings


class StartSamplingExperimentFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self, experiment_id: str):
        configuration = Configuration(
            host=self._settings.reinvent_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            api_instance.sampling_jobs_job_id_sampling_post(job_id=experiment_id)
