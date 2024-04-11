import msa_light_microservice as microservice

from msa_light_microservice import DefaultApi, ApiClient, Configuration
from nolabs.api_models.drug_discovery import PredictMsaRequest, PredictMsaResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings


class GenerateMsaFeature:
    def __init__(self, file_management: TargetsFileManagement,
                 settings: Settings):
        self._file_management = file_management
        self._settings = settings

    def handle(self, request: PredictMsaRequest) -> PredictMsaResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        fasta_contents = self._file_management.get_fasta_contents(experiment_id, target_id)

        configuration = Configuration(
            host=self._settings.msa_light_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            request = microservice.RunMsaPredictionRequest(
                api_url=self._settings.drug_discovery_msa_remote_prediction_url,
                fasta_contents=fasta_contents)
            msa_contents = api_instance.predict_msa_predict_msa_post(predict_msa_predict_msa_post=request).msa_contents

        self._file_management.store_msa(experiment_id, target_id, msa_contents)

        return PredictMsaResponse(msa_contents)
