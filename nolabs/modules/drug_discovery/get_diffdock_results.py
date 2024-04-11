from nolabs.api_models.drug_discovery import GetDiffDockDockingResultDataResponse, GetDockingResultDataRequest, \
    DiffDockLigandMetaData, GetDiffDockLigandSdfRequest, GetDiffDockLigandSdfResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.result import JobId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement


class GetDiffDockDockingResultsFeature:
    def __init__(self, result_file_management: ResultsFileManagement):
        self._result_file_management = result_file_management

    def handle(self, request: GetDockingResultDataRequest) -> GetDiffDockDockingResultDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        predicted_pdb, result_ligands_list = self._result_file_management.get_diffdock_docking_result_data(
            experiment_id, target_id, ligand_id, job_id)
        predicted_ligands = [DiffDockLigandMetaData(
            job_id=job_id.value,
            ligand_id=ligand_id.value,
            target_id=target_id.value,
            predicted_ligand_file_name=result_ligand.ligand_file_name,
            confidence=result_ligand.confidence,
            minimized_affinity=result_ligand.minimized_affinity,
            scored_affinity=result_ligand.scored_affinity
        ) for result_ligand in result_ligands_list]

        return GetDiffDockDockingResultDataResponse(predicted_pdb=predicted_pdb, predicted_ligands=predicted_ligands)


class GetDiffDockLigandSdfFeature:
    def __init__(self, result_file_management: ResultsFileManagement):
        self._result_file_management = result_file_management

    def handle(self, request: GetDiffDockLigandSdfRequest) -> GetDiffDockLigandSdfResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)
        ligand_file_name = request.ligand_file_name

        ligand_sdf = self._result_file_management.get_diffdock_ligand_data(experiment_id,
                                                                           target_id,
                                                                           ligand_id,
                                                                           job_id,
                                                                           ligand_file_name).predicted_sdf_contents

        return GetDiffDockLigandSdfResponse(sdf_contents=ligand_sdf)
