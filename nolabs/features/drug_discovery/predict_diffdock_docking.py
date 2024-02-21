from typing import List

from nolabs.features.drug_discovery.data_models.ligand import LigandId
from nolabs.features.drug_discovery.data_models.result import JobId, UmolDockingResultData, DiffDockDockingResultData
from nolabs.infrastructure.settings import Settings

import diffdock_microservice
from diffdock_microservice import DefaultApi as DiffDockDefaultApi
from diffdock_microservice import ApiClient as DiffDockApiClient
from diffdock_microservice import Configuration as DiffDockConfiguration

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

        pdb_content = self.get_folding(experiment_id, target_id, job_id, job_metadata.folding_method)
        _, sdf_contents, _ = self._ligand_file_management.get_ligand_data(experiment_id, target_id, ligand_id)

        configuration = DiffDockConfiguration(
            host=self._settings.diffdock_host,
        )
        with DiffDockApiClient(configuration=configuration) as client:
            api_instance = DiffDockDefaultApi(client)
            request = diffdock_microservice.RunDiffDockPredictionRequest(pdb_contents=pdb_content,
                                                                         sdf_contents=sdf_contents,
                                                                         job_id=job_id.value,
                                                                         )
            response = api_instance.predict_run_docking_post(run_diff_dock_prediction_request=request)

        predicted_pdb = response.pdb_contents
        predicted_ligands_list = []
        for ligand_result in response.sdf_results:
            confidence = ligand_result.confidence
            minimised_affinity = ligand_result.minimized_affinity
            scored_affinity = ligand_result.scored_affinity
            ligand_file_name = ligand_result.ligand_file_name
            ligand_sdf = ligand_result.ligand_sdf
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
                                                                 ligand_file_name=result_data.predicted_sdf_file_name))

        return RunDiffDockDockingJobResponse(predicted_pdb=predicted_pdb,
                                             predicted_ligands=predicted_ligands_list)

    def get_folding(self, experiment_id: ExperimentId, target_id: TargetId, job_id: JobId, folding_method: str) -> str:
        if not self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method):
            _, sequence, _ = self._target_file_management.get_target_data(experiment_id, target_id)

            host = None

            if folding_method == "esmfold_light":
                host = self._settings.esmfold_light_host
            else:
                host = self._settings.esmfold_host

            configuration = FoldingConfiguration(
                host=host,
            )
            with FoldingApiClient(configuration=configuration) as client:
                api_instance = FoldingDefaultApi(client)
                request = folding_microservice.RunEsmFoldPredictionRequest(protein_sequence=sequence,
                                                                           job_id=job_id.value)
                pdb_contents = api_instance.predict_through_api_run_folding_post(
                    run_esm_fold_prediction_request=request).pdb_content
            self._target_file_management.store_pdb_contents(experiment_id, target_id, pdb_contents, folding_method)
            return pdb_contents
        else:
            return self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method)
