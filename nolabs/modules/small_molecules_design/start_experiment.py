from reinvent_microservice import Configuration, ApiClient, DefaultApi

from nolabs.infrastructure.settings import Settings


class StartExperimentFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self, experiment_id: str):
        configuration = Configuration(
            host=self._settings.reinvent_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            api_instance.run_jobs_job_id_run_post(job_id=experiment_id)
