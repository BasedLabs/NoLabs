import diffdock_microservice
from diffdock_microservice import ApiClient as DiffDockApiClient
from diffdock_microservice import Configuration as DiffDockConfiguration
from diffdock_microservice import DefaultApi as DiffDockDefaultApi

from nolabs.api_models.drug_discovery import RunDiffDockDockingJobRequest, \
    RunDiffDockDockingJobResponse, DiffDockLigandMetaData
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.result import JobId, DiffDockDockingResultData
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.folding_backends_factory import run_folding
from nolabs.modules.drug_discovery.services.folding_methods import FoldingMethods
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings


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
        _, sdf_contents, _ = self._ligand_file_management.get_target_ligand_data(experiment_id, target_id, ligand_id)

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

    def get_folding(self, experiment_id: ExperimentId, target_id: TargetId, folding_method: FoldingMethods) -> str:
        if not self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method):
            _, sequence, _ = self._target_file_management.get_target_data(experiment_id, target_id)

            FoldingMethods.ensure_contains(folding_method)

            pdb_content = run_folding(sequence, self._settings, folding_method)
            self._target_file_management.store_pdb_contents(experiment_id, target_id, pdb_content, folding_method)
            self._target_file_management.update_target_metadata(experiment_id, target_id, "folding_method",
                                                                folding_method)

        return self._target_file_management.get_pdb_contents(experiment_id, target_id, folding_method)