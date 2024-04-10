import datetime
from typing import List

from reinvent_microservice import Configuration, ApiClient, ReinventApi

from nolabs.api_models.small_molecules_design import SmilesResponse
from nolabs.infrastructure.settings import Settings
from nolabs.modules.small_molecules_design.services.file_management import FileManagement


class GetSmilesFeature:
    def __init__(self, file_management: FileManagement, settings: Settings):
        self._settings = settings
        self._fm = file_management

    async def handle(self, experiment_id: str) -> List[SmilesResponse]:
        if not self._fm.experiment_exists(experiment_id):
            return []

        configuration = Configuration(
            host=self._settings.reinvent_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = ReinventApi(client)

            response = api_instance.smiles_api_reinvent_config_id_smiles_get(experiment_id)

            return [SmilesResponse(
                smiles=s.smiles,
                drug_likeness=s.drug_likeness,
                score=s.score,
                stage=s.stage,
                created_at=datetime.datetime.utcnow()
            ) for s in response.smiles]
