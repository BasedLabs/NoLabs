from reinvent_microservice import Configuration, ApiClient, DefaultApi

from nolabs.domain.experiment import ExperimentId
from nolabs.modules.small_molecules_design.services.file_management import FileManagement
from nolabs.infrastructure.settings import Settings


class DeleteExperimentFeature:
    def __init__(self, settings: Settings, file_management: FileManagement):
        self._settings = settings
        self._fm = file_management

    async def handle(self, experiment_id: str):
        self._fm.delete_experiment_folder(ExperimentId(experiment_id))
        configuration = Configuration(
            host=self._settings.reinvent_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            api_instance.delete_jobs_job_id_delete(experiment_id)

