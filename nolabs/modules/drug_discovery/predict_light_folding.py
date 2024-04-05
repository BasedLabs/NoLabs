from nolabs.modules.drug_discovery.services.folding_backends_factory import run_folding
from nolabs.api_models.drug_discovery import PredictFoldingRequest, PredictFoldingResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.folding_methods import FoldingMethods
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings


class PredictEsmFoldLightFeature:
    def __init__(self, file_management: TargetsFileManagement,
                 settings: Settings):
        self._file_management = file_management
        self._settings = settings

    def handle(self, request: PredictFoldingRequest) -> PredictFoldingResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        _, sequence, _ = self._file_management.get_target_data(experiment_id, target_id)

        if len(sequence) > 400:
            raise NoLabsException(messages=["Light folding does not support sequences longer than 400. Please use "
                                            "other folding backend"],
                                  error_code=ErrorCodes.drug_discovery_folding_error)

        pdb_content = run_folding(sequence, self._settings, FoldingMethods.esmfold_light)

        self._file_management.store_pdb_contents(experiment_id, target_id, pdb_content, FoldingMethods.esmfold_light)
        self._file_management.update_target_metadata(experiment_id, target_id, "folding_method", FoldingMethods.esmfold_light)

        return PredictFoldingResponse(pdb_content=pdb_content)
