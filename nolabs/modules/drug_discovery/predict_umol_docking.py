from typing import List

# TODO: Add self-hosted MSA microservice
# if settings.drug_discovery_self_hosted_msa:
import msa_light_microservice as msa_microservice
import p2rank_microservice
import umol_microservice
from msa_light_microservice import ApiClient as MsaApiClient
from msa_light_microservice import Configuration as MsaConfiguration
from msa_light_microservice import DefaultApi as MsaDefaultApi
from p2rank_microservice import ApiClient as P2RankApiClient
from p2rank_microservice import Configuration as P2RankConfiguration
from p2rank_microservice import DefaultApi as P2RankDefaultApi
from umol_microservice import ApiClient as UmolApiClient
from umol_microservice import Configuration as UmolConfiguration
from umol_microservice import DefaultApi as UmolDefaultApi

from nolabs.api_models.drug_discovery import RunUmolDockingJobRequest, RunUmolDockingJobResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.result import JobId, UmolDockingResultData
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.folding_backends_factory import run_folding
from nolabs.modules.drug_discovery.services.folding_methods import FoldingMethods
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings


class PredictUmolDockingFeature:
    def __init__(self, target_file_management: TargetsFileManagement,
                 ligand_file_management: LigandsFileManagement,
                 result_file_management: ResultsFileManagement,
                 settings: Settings):
        self._target_file_management = target_file_management
        self._ligand_file_management = ligand_file_management
        self._result_file_management = result_file_management
        self._settings = settings

    def handle(self, request: RunUmolDockingJobRequest) -> RunUmolDockingJobResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        job_metadata = self._result_file_management.get_job_metadata(experiment_id, target_id, ligand_id, job_id)

        msa_contents = self.get_msa(experiment_id, target_id, job_id)
        pocket_ids = self.get_pocket_ids(experiment_id, target_id, job_id, job_metadata.folding_method)
        self._result_file_management.store_result_input_pocketIds(experiment_id,
                                                                  target_id,
                                                                  ligand_id,
                                                                  job_id,
                                                                  pocket_ids[:])
        fasta_contents = self._target_file_management.get_fasta_contents(experiment_id, target_id)
        _, _, ligand_smiles = self._ligand_file_management.get_ligand_data(experiment_id, target_id, ligand_id)

        configuration = UmolConfiguration(
            host=self._settings.umol_host,
        )
        with UmolApiClient(configuration=configuration) as client:
            api_instance = UmolDefaultApi(client)
            request = umol_microservice.RunUmolPredictionRequest(protein_sequence=fasta_contents,
                                                                 msa_content=msa_contents,
                                                                 pocket_ids=pocket_ids,
                                                                 ligand_smiles=ligand_smiles,
                                                                 job_id=job_id.value)
            response = api_instance.predict_run_umol_post(run_umol_prediction_request=request)

        predicted_pdb = response.pdb_contents
        predicted_sdf = response.sdf_contents
        plddt_array = response.plddt_array

        result_data = UmolDockingResultData(predicted_pdb=predicted_pdb,
                                            predicted_sdf=predicted_sdf,
                                            plddt_array=plddt_array)

        self._result_file_management.store_umol_result_data(experiment_id, target_id, ligand_id, job_id, result_data)

        return RunUmolDockingJobResponse(predicted_pdb=predicted_pdb,
                                         predicted_sdf=predicted_sdf,
                                         plddt_array=plddt_array,
                                         job_id=job_id.value)

    def get_folding(self, experiment_id: ExperimentId, target_id: TargetId, folding_method: FoldingMethods) -> str:
        if not self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method):
            _, sequence, _ = self._target_file_management.get_target_data(experiment_id, target_id)

            FoldingMethods.ensure_contains(folding_method)

            pdb_content = run_folding(sequence, self._settings, folding_method)
            self._target_file_management.store_pdb_contents(experiment_id, target_id, pdb_content, folding_method)
            self._target_file_management.update_target_metadata(experiment_id, target_id, "folding_method",
                                                                folding_method)

        return self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method)

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
                    job_id=job_id.value
                )
                msa_contents = api_instance.predict_msa_predict_msa_post(
                    run_msa_prediction_request=request).msa_contents
            self._target_file_management.store_msa(experiment_id, target_id, msa_contents)
            return msa_contents
        else:
            return self._target_file_management.get_msa(experiment_id, target_id)

    def get_pocket_ids(self, experiment_id: ExperimentId,
                       target_id: TargetId,
                       job_id: JobId,
                       folding_method: str) -> List[int]:
        if not self._target_file_management.get_binding_pocket(experiment_id, target_id):
            configuration = P2RankConfiguration(
                host=self._settings.p2rank_host,
            )
            with P2RankApiClient(configuration=configuration) as client:
                api_instance = P2RankDefaultApi(client)
                pdb_contents = self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method)
                pdb_fixer_request = p2rank_microservice.RunP2RankPredictionRequest(pdb_contents=pdb_contents,
                                                                                   job_id=job_id.value)
                pocket_ids = api_instance.predict_run_p2rank_post(
                    run_p2_rank_prediction_request=pdb_fixer_request).pocket_ids
                self._target_file_management.store_binding_pocket(experiment_id, target_id, pocket_ids)
            return pocket_ids
        else:
            return self._target_file_management.get_binding_pocket(experiment_id, target_id)
