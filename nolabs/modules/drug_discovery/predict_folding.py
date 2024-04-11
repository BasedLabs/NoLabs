from nolabs.api_models.drug_discovery import PredictFoldingRequest, PredictFoldingResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.folding_backends_factory import run_folding
from nolabs.modules.drug_discovery.services.folding_methods import FoldingMethods
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings


class PredictEsmFoldFeature:
    def __init__(self, file_management: TargetsFileManagement,
                 settings: Settings):
        self._file_management = file_management
        self._settings = settings

    def handle(self, request: PredictFoldingRequest) -> PredictFoldingResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        _, sequence, _ = self._file_management.get_target_data(experiment_id, target_id)

        pdb_content = run_folding(sequence, self._settings, FoldingMethods.esmfold)

        self._file_management.store_pdb_contents(experiment_id, target_id, pdb_content, FoldingMethods.esmfold)
        self._file_management.update_target_metadata(experiment_id, target_id, "folding_method", FoldingMethods.esmfold)

        return PredictFoldingResponse(pdb_content=pdb_content)
