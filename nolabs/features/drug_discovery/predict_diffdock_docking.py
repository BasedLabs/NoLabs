from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.features.drug_discovery.data_models.ligand import LigandId
from nolabs.features.drug_discovery.data_models.result import JobId, DiffDockDockingResultData
from nolabs.infrastructure.settings import Settings

import diffdock_microservice
from diffdock_microservice import DefaultApi as DiffDockDefaultApi
from diffdock_microservice import ApiClient as DiffDockApiClient
from diffdock_microservice import Configuration as DiffDockConfiguration

import esmfold_light_microservice
from esmfold_light_microservice import DefaultApi as EsmFoldLightDefaultApi
from esmfold_light_microservice import ApiClient as EsmFoldLightApiClient
from esmfold_light_microservice import Configuration as EsmFoldLightConfiguration

import esmfold_microservice
from esmfold_microservice import DefaultApi as EsmFoldDefaultApi
from esmfold_microservice import ApiClient as EsmFoldApiClient
from esmfold_microservice import Configuration as EsmFoldConfiguration

from nolabs.api_models.drug_discovery import RunDiffDockDockingJobRequest, \
    RunDiffDockDockingJobResponse, DiffDockLigandMetaData
from nolabs.domain.experiment import ExperimentId
from nolabs.features.drug_discovery.data_models.target import TargetId
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.features.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.features.drug_discovery.services.result_file_management import ResultsFileManagement


class PredictDiffDockDockingFeature:
    def __init__(self, target_file_management: TargetsFileManagement,
                 ligand_file_management: LigandsFileManagement,
                 result_file_management: ResultsFileManagement,
                 settings: Settings):
        self._target_file_management = target_file_management
        self._ligand_file_management = ligand_file_management
        self._result_file_management = result_file_management
        self._settings = settings

    def handle(self, request: RunDiffDockDockingJobRequest) -> RunDiffDockDockingJobResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        job_metadata = self._result_file_management.get_job_metadata(experiment_id, target_id, ligand_id, job_id)
        params = self._result_file_management.get_diffdock_params(experiment_id, target_id, ligand_id, job_id)

        pdb_content = self.get_folding(experiment_id, target_id, job_metadata.folding_method)
        _, sdf_contents, _ = self._ligand_file_management.get_ligand_data(experiment_id, target_id, ligand_id)

        configuration = DiffDockConfiguration(
            host=self._settings.diffdock_host,
        )
        with DiffDockApiClient(configuration=configuration) as client:
            api_instance = DiffDockDefaultApi(client)
            request = diffdock_microservice.RunDiffDockPredictionRequest(pdb_contents=pdb_content,
                                                                         sdf_contents=sdf_contents,
                                                                         samples_per_complex=params.samples_per_complex,
                                                                         job_id=job_id.value,
                                                                         )
            response = api_instance.predict_run_docking_post(run_diff_dock_prediction_request=request)

        predicted_pdb = response.pdb_contents
        predicted_ligands_list = []
        for ligand_result in response.sdf_results:
            confidence = ligand_result.confidence
            if confidence:
                minimised_affinity = ligand_result.minimized_affinity
                scored_affinity = ligand_result.scored_affinity
                ligand_file_name = ligand_result.sdf_file_name
                ligand_sdf = ligand_result.sdf_content
                result_data = DiffDockDockingResultData(confidence=confidence,
                                                        scored_affinity=scored_affinity,
                                                        minimized_affinity=minimised_affinity,
                                                        predicted_sdf_file_name=ligand_file_name,
                                                        predicted_sdf_contents=ligand_sdf,
                                                        predicted_pdb=predicted_pdb)

                self._result_file_management.store_diffdock_result_data(experiment_id,
                                                                        target_id,
                                                                        ligand_id,
                                                                        job_id,
                                                                        result_data)

                predicted_ligands_list.append(DiffDockLigandMetaData(target_id=target_id.value,
                                                                     ligand_id=ligand_id.value,
                                                                     job_id=job_id.value,
                                                                     scored_affinity=result_data.scored_affinity,
                                                                     minimized_affinity=result_data.minimized_affinity,
                                                                     confidence=result_data.confidence,
                                                                     predicted_ligand_file_name=result_data.predicted_sdf_file_name))

        return RunDiffDockDockingJobResponse(predicted_pdb=predicted_pdb,
                                             predicted_ligands=predicted_ligands_list)

    def get_folding(self, experiment_id: ExperimentId, target_id: TargetId, folding_method: str) -> str:
        if not self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method):
            _, sequence, _ = self._target_file_management.get_target_data(experiment_id, target_id)

            pdb_contents = None

            if folding_method == "esmfold_light":
                pdb_contents = self.predict_esmfold_light(experiment_id, target_id)
            else:
                pdb_contents = self.predict_esmfold(experiment_id, target_id)

            return pdb_contents
        else:
            return self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method)

    def predict_esmfold_light(self, experiment_id: ExperimentId, target_id: TargetId):

        configuration = EsmFoldLightConfiguration(
            host=self._settings.esmfold_light_host,
        )
        with EsmFoldLightApiClient(configuration=configuration) as client:
            api_instance = EsmFoldLightDefaultApi(client)
            _, sequence, _ = self._target_file_management.get_target_data(experiment_id, target_id)
            if len(sequence) > 400:
                raise NoLabsException(messages=["Light folding does not support sequences longer than 400. Please use "
                                                "other folding backend"],
                                      error_code=ErrorCodes.drug_discovery_folding_error)
            else:
                request = esmfold_light_microservice.RunEsmFoldPredictionRequest(protein_sequence=sequence)
                pdb_content = api_instance.predict_through_api_run_folding_post(
                    run_esm_fold_prediction_request=request).pdb_content

                self._target_file_management.store_pdb_contents(experiment_id, target_id, pdb_content, "esmfold_light")
                self._target_file_management.update_target_metadata(experiment_id, target_id, "folding_method",
                                                             "esmfold_light")
                return pdb_content


    def predict_esmfold(self, experiment_id: ExperimentId, target_id: TargetId):

        configuration = EsmFoldConfiguration(
            host=self._settings.esmfold_host,
        )
        with EsmFoldApiClient(configuration=configuration) as client:
            api_instance = EsmFoldDefaultApi(client)
            _, sequence, _ = self._target_file_management.get_target_data(experiment_id, target_id)

            request = esmfold_microservice.RunEsmFoldPredictionRequest(protein_sequence=sequence)
            pdb_content = api_instance.predict_run_folding_post(
                run_esm_fold_prediction_request=request).pdb_content

            self._target_file_management.store_pdb_contents(experiment_id, target_id, pdb_content, "esmfold")
            self._target_file_management.update_target_metadata(experiment_id, target_id, "folding_method",
                                                         "esmfold")
            return pdb_content