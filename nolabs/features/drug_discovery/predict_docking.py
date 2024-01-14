from typing import List

from nolabs.infrastructure.settings import Settings
import p2rank_microservice
import umol_microservice

settings = Settings()
if settings.is_light_infrastructure:
    import esmfold_light_microservice as folding_microservice
else:
    import esmfold_microservice as folding_microservice

# TODO: Add self-hosted MSA microservice
# if settings.drug_discovery_self_hosted_msa:
import msa_light_microservice as msa_microservice

from nolabs.api_models.drug_discovery import DockingRequest, DockingResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.features.drug_discovery.services.ligand_file_management import LigandsFileManagement


class PredictDockingFeature:
    def __init__(self, target_file_management: TargetsFileManagement, ligand_file_management: LigandsFileManagement):
        self._target_file_management = target_file_management
        self._ligand_file_management = ligand_file_management

    def handle(self, request: DockingRequest) -> DockingResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        pdb_contents = self.get_folding(experiment_id, target_id)
        msa_contents = self.get_msa(experiment_id, target_id)
        pocket_ids = self.get_pocket_ids(experiment_id, target_id)

        return DockingResponse(predicted_pdb=pdb_contents, predicted_sdf="")

    def get_folding(self, experiment_id: ExperimentId, target_id: TargetId) -> str:
        if not self._target_file_management.get_pdb_contents(experiment_id, target_id):
            client = folding_microservice.DefaultApi(api_client=folding_microservice.ApiClient())
            _, sequence, _ = self._target_file_management.get_target_data(experiment_id, target_id)
            request = folding_microservice.RunEsmFoldPredictionRequest(protein_sequence=sequence)
            pdb_contents = client.predict_through_api_run_folding_post(request=request).pdb_content
            self._target_file_management.store_pdb_contents(experiment_id, target_id, pdb_contents)
            return pdb_contents
        else:
            return self._target_file_management.get_pdb_contents(experiment_id, target_id)

    def get_msa(self, experiment_id: ExperimentId, target_id: TargetId) -> str:
        if not self._target_file_management.get_msa(experiment_id, target_id):
            client = msa_microservice.DefaultApi(api_client=msa_microservice.ApiClient())
            fasta_contents = self._target_file_management.get_fasta_contents(experiment_id, target_id)
            msa_server_url = self._target_file_management.get_msa_api_url()
            request = msa_microservice.RunMsaPredictionRequest(api_url=msa_server_url, fasta_contents=fasta_contents)
            msa_contents = client.predict_msa_predict_msa_post(request=request).msa_contents
            self._target_file_management.store_msa(experiment_id, target_id, msa_contents)
            return msa_contents
        else:
            return self._target_file_management.get_msa(experiment_id, target_id)

    def get_pocket_ids(self, experiment_id: ExperimentId, target_id: TargetId) -> List[int]:
        if not self._target_file_management.get_binding_pocket(experiment_id, target_id):
            client = p2rank_microservice.DefaultApi(api_client=p2rank_microservice.ApiClient())
            pdb_contents = self._target_file_management.get_pdb_contents(experiment_id, target_id)
            pdb_fixer_request = p2rank_microservice.RunP2RankPredictionRequest(pdb_contents=pdb_contents)
            pocket_ids = client.predict_run_p2rank_post(request=pdb_fixer_request).pocket_ids
            self._target_file_management.store_binding_pocket(experiment_id, target_id, pocket_ids)
            return pocket_ids
        else:
            return self._target_file_management.get_binding_pocket(experiment_id, target_id)
