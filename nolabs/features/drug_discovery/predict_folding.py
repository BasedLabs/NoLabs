import esmfold_microservice as microservice

from esmfold_microservice import ApiClient, DefaultApi, Configuration

from nolabs.api_models.drug_discovery import PredictFoldingRequest, PredictFoldingResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement
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

        configuration = Configuration(
            host=self._settings.esmfold_host,
        )
        print(self._settings.esmfold_host)
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            _, sequence, _ = self._file_management.get_target_data(experiment_id, target_id)
            request = microservice.RunEsmFoldPredictionRequest(protein_sequence=sequence)
            pdb_content = api_instance.predict_run_folding_post(run_esm_fold_prediction_request=request).pdb_content

            self._file_management.store_pdb_contents(experiment_id, target_id, pdb_content)

        return PredictFoldingResponse(pdb_content=pdb_content)
