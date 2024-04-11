from nolabs.api_models.drug_discovery import GetDockingParamsRequest, GetDockingParamsResponse, \
    UpdateDockingParamsRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.result import JobId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement


class GetDockingParamsFeature:
    def __init__(self, file_management: ResultsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetDockingParamsRequest) -> GetDockingParamsResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        job_metadata = self._file_management.get_job_metadata(experiment_id, target_id, ligand_id, job_id)

        return GetDockingParamsResponse(folding_method=job_metadata.folding_method, docking_method=job_metadata.docking_method)


class UpdateDockingParamsFeature:
    def __init__(self, file_management: ResultsFileManagement):
        self._file_management = file_management

    def handle(self, request: UpdateDockingParamsRequest) -> GetDockingParamsResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)
        folding_method = request.folfing_method
        docking_method = request.docking_method

        self._file_management.update_job_metadata(experiment_id, target_id, ligand_id, job_id, "folding_method", folding_method)
        self._file_management.update_job_metadata(experiment_id, target_id, ligand_id, job_id, "docking_method",
                                                  docking_method)

        return GetDockingParamsResponse(folding_method=folding_method, docking_method=docking_method)
