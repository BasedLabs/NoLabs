import msa_light_microservice as microservice

from nolabs.api_models.drug_discovery import PredictMsaRequest, PredictMsaResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement

class GenerateMsaFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: PredictMsaRequest) -> PredictMsaResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        client = microservice.DefaultApi(api_client=microservice.ApiClient())
        fasta_contents = self._file_management.get_fasta_contents(experiment_id, target_id)
        msa_server_url = self._file_management.get_msa_api_url()
        request = microservice.RunMsaPredictionRequest(api_url=msa_server_url, fasta_contents=fasta_contents)
        msa_contents = client.predict_msa_predict_msa_post(request=request).msa_contents

        self._file_management.store_msa(experiment_id, target_id, msa_contents)

        return PredictMsaResponse(msa_contents)