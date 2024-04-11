from nolabs.api_models.drug_discovery import DeleteDockingJobRequest, DeleteDockingJobResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.result import JobId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement

class DeleteJobFeature:
    def __init__(self, file_management: ResultsFileManagement):
        self._file_management = file_management

    def handle(self, request: DeleteDockingJobRequest) -> DeleteDockingJobResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        deleted_job_id = self._file_management.delete_result(experiment_id, target_id, ligand_id, job_id)

        return DeleteDockingJobResponse(job_id=deleted_job_id.value)