from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.result import JobId

from nolabs.api_models.drug_discovery import GetDockingResultDataRequest, GetUmolDockingResultDataResponse, \
    GetDiffDockDockingResultDataResponse, DiffDockLigandMetaData, GetDiffDockLigandSdfRequest, \
    GetDiffDockLigandSdfResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement


class GetUmolDockingResultsFeature:
    def __init__(self, result_file_management: ResultsFileManagement):
        self._result_file_management = result_file_management

    def handle(self, request: GetDockingResultDataRequest) -> GetUmolDockingResultDataResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        result_data = self._result_file_management.get_umol_docking_result_data(experiment_id, target_id, ligand_id,
                                                                                job_id)
        return GetUmolDockingResultDataResponse(predicted_pdb=result_data.predicted_pdb,
                                                predicted_sdf=result_data.predicted_sdf,
                                                plddt_array=result_data.plddt_array,
                                                job_id=job_id.value)


