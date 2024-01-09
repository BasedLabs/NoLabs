import conformations_microservice as microservice
from conformations_microservice import RunPdbSimulationsRequest, ForceFields, WaterForceFields

import nolabs.api_models.protein_design as api
from nolabs.infrastructure.settings import Settings
from nolabs.features.protein_design.services.file_management import FileManagement
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.infrastructure.websockets import ConformationsWebsocket
from nolabs.exceptions import NoLabsException, ErrorCodes


class RunProteinDesignFeature:
    async def handle(self,
                     request: api.RunProteinDesignRequest) -> api.RunProteinDesignResponse:
        assert request

        ex
