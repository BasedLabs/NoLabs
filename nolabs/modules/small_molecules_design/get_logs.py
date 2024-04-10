from typing import List

from reinvent_microservice import Configuration, ApiClient, ReinventApi

from nolabs.api_models.small_molecules_design import LogsResponse
from nolabs.modules.small_molecules_design.services.file_management import FileManagement
from nolabs.infrastructure.settings import Settings


class GetLogsFeature:
    def __init__(self, file_management: FileManagement, settings: Settings):
        self._settings = settings
        self._fm = file_management

    async def handle(self, experiment_id: str) -> LogsResponse:
        configuration = Configuration(
            host=self._settings.reinvent_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = ReinventApi(client)

            response = api_instance.logs_api_reinvent_config_id_logs_get(experiment_id)

            if not response:
                return LogsResponse(
                    output='Nothing',
                    docking_output='Nothing',
                    errors='Nothing'
                )

            instance = response.actual_instance

            return LogsResponse(output=instance.output,
                                docking_output=instance.docking_output,
                                errors=instance.errors)
