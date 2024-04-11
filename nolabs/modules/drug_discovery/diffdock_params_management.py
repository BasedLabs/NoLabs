from nolabs.api_models.drug_discovery import GetDiffDockParamsRequest, GetDiffDockParamsResponse, \
    UpdateDiffDockParamsRequest
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.result import JobId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement


class GetDiffDockParamsFeature:
    def __init__(self, file_management: ResultsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetDiffDockParamsRequest) -> GetDiffDockParamsResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        params = self._file_management.get_diffdock_params(experiment_id, target_id, ligand_id, job_id)

        return GetDiffDockParamsResponse(samples_per_complex=params.samples_per_complex)


class UpdateDiffDockParamsFeature:
    def __init__(self, file_management: ResultsFileManagement):
        self._file_management = file_management

    def handle(self, request: UpdateDiffDockParamsRequest) -> GetDiffDockParamsResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)
        new_samples_per_complex = request.samples_per_complex

        params = self._file_management.get_diffdock_params(experiment_id, target_id, ligand_id, job_id)
        params.samples_per_complex = new_samples_per_complex

        self._file_management.update_diffdock_params(experiment_id, target_id, ligand_id, job_id, params)

        return GetDiffDockParamsResponse(samples_per_complex=params.samples_per_complex)
