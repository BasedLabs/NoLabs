from typing import List

from nolabs.features.drug_discovery.data_models.result import JobId
from nolabs.infrastructure.settings import Settings

# TODO: Add self-hosted MSA microservice
# if settings.drug_discovery_self_hosted_msa:
import msa_light_microservice as msa_microservice
from msa_light_microservice import DefaultApi as MsaDefaultApi
from msa_light_microservice import ApiClient as MsaApiClient
from msa_light_microservice import Configuration as MsaConfiguration

import p2rank_microservice
from p2rank_microservice import DefaultApi as P2RankDefaultApi
from p2rank_microservice import ApiClient as P2RankApiClient
from p2rank_microservice import Configuration as P2RankConfiguration

import umol_microservice
from umol_microservice import DefaultApi as UmolDefaultApi
from umol_microservice import ApiClient as UmolApiClient
from umol_microservice import Configuration as UmolConfiguration

settings = Settings()
if settings.is_light_infrastructure:
    import esmfold_light_microservice as folding_microservice
    from esmfold_light_microservice import DefaultApi as FoldingDefaultApi
    from esmfold_light_microservice import ApiClient as FoldingApiClient
    from esmfold_light_microservice import Configuration as FoldingConfiguration
else:
    import esmfold_microservice as folding_microservice
    from esmfold_microservice import DefaultApi as FoldingDefaultApi
    from esmfold_microservice import ApiClient as FoldingApiClient
    from esmfold_microservice import Configuration as FoldingConfiguration

from nolabs.api_models.drug_discovery import CheckJobIsRunningRequest, CheckJobIsRunningResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.target import TargetId

class CheckUmolRunningFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self, request: CheckJobIsRunningRequest) -> CheckJobIsRunningResponse:
        assert request

        job_id = JobId(request.job_id)

        configuration = UmolConfiguration(
            host=self._settings.umol_host,
        )
        with UmolApiClient(configuration=configuration) as client:
            api_instance = UmolDefaultApi(client)
            response = api_instance.is_job_running_job_job_id_is_running_get(job_id=job_id.value)

        return CheckJobIsRunningResponse(is_running=False)

class CheckMsaRunningFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self, request: CheckJobIsRunningRequest) -> CheckJobIsRunningResponse:

        job_id = JobId(request.job_id)

        configuration = MsaConfiguration(
            host=self._settings.msa_light_host,
        )
        with MsaApiClient(configuration=configuration) as client:
            api_instance = MsaDefaultApi(client)
            response = api_instance.is_job_running_job_job_id_is_running_get(job_id=job_id.value).is_running

        return CheckJobIsRunningResponse(is_running=response)

    def get_folding(self, experiment_id: ExperimentId, target_id: TargetId, job_id: JobId) -> str:
        if not self._target_file_management.get_pdb_contents(experiment_id, target_id):
            _, sequence, _ = self._target_file_management.get_target_data(experiment_id, target_id)
            configuration = FoldingConfiguration(
                host=self._settings.folding_host,
            )
            with FoldingApiClient(configuration=configuration) as client:
                api_instance = FoldingDefaultApi(client)
                request = folding_microservice.RunEsmFoldPredictionRequest(protein_sequence=sequence)
                pdb_contents = api_instance.predict_through_api_run_folding_post(run_esm_fold_prediction_request=request).pdb_content
            self._target_file_management.store_pdb_contents(experiment_id, target_id, pdb_contents)
            return pdb_contents
        else:
            return self._target_file_management.get_pdb_contents(experiment_id, target_id)

    def get_msa(self, experiment_id: ExperimentId, target_id: TargetId, job_id: JobId) -> str:
        if not self._target_file_management.get_msa(experiment_id, target_id):

            fasta_contents = self._target_file_management.get_fasta_contents(experiment_id, target_id)

            configuration = MsaConfiguration(
                host=self._settings.msa_light_host,
            )
            with MsaApiClient(configuration=configuration) as client:
                api_instance = MsaDefaultApi(client)
                request = msa_microservice.RunMsaPredictionRequest(
                    api_url=self._settings.drug_discovery_msa_remote_prediction_url,
                    fasta_contents=fasta_contents,
                )
                msa_contents = api_instance.predict_msa_predict_msa_post(run_msa_prediction_request=request).msa_contents
            self._target_file_management.store_msa(experiment_id, target_id, msa_contents)
            return msa_contents
        else:
            return self._target_file_management.get_msa(experiment_id, target_id)

    def get_pocket_ids(self, experiment_id: ExperimentId, target_id: TargetId, job_id: JobId) -> List[int]:
        if not self._target_file_management.get_binding_pocket(experiment_id, target_id):
            configuration = P2RankConfiguration(
                host=self._settings.p2rank_host,
            )
            with P2RankApiClient(configuration=configuration) as client:
                api_instance = P2RankDefaultApi(client)
                pdb_contents = self._target_file_management.get_pdb_contents(experiment_id, target_id)
                pdb_fixer_request = p2rank_microservice.RunP2RankPredictionRequest(pdb_contents=pdb_contents)
                pocket_ids = api_instance.predict_run_p2rank_post(run_p2_rank_prediction_request=pdb_fixer_request).pocket_ids
                self._target_file_management.store_binding_pocket(experiment_id, target_id, pocket_ids)
            return pocket_ids
        else:
            return self._target_file_management.get_binding_pocket(experiment_id, target_id)
