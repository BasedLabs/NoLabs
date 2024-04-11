from nolabs.modules.drug_discovery.data_models.result import JobId
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.modules.drug_discovery.services.result_file_management import ResultsFileManagement
from nolabs.api_models.drug_discovery import GetResultsListForTargetLigandRequest, \
    GetResultsListForTargetLigandResponse, \
    GetAllResultsListRequest, GetAllResultsListResponse, ResultMetaData, \
    CheckResultDataAvailableRequest, CheckResultDataAvailableResponse, \
    CheckMsaDataAvailableResponse, CheckMsaDataAvailableRequest, \
    CheckPocketDataAvailableRequest, CheckPocketDataAvailableResponse, CheckFoldingDataAvailableResponse, \
    CheckFoldingDataAvailableRequest, GetJobBindingPocketDataRequest, GetJobBindingPocketDataResponse, \
    GetAllJobsListRequest, GetAllJobsListResponse, JobMetaData, GetJobsListForTargetLigandRequest, \
    GetJobsListForTargetLigandResponse
from nolabs.domain.experiment import ExperimentId
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.services.target_file_management import TargetsFileManagement


class GetAllResultsListFeature:
    def __init__(self, results_file_management: ResultsFileManagement,
                 target_file_management: TargetsFileManagement,
                 ligand_file_management: LigandsFileManagement):
        self._results_file_management = results_file_management
        self._targets_file_management = target_file_management
        self._ligands_file_management = ligand_file_management

    def handle(self, request: GetAllResultsListRequest) -> GetAllResultsListResponse:
        experiment_id = ExperimentId(request.experiment_id)

        results_list = []

        for target in self._targets_file_management.get_targets_list(experiment_id):
            target_id = TargetId(target.target_id)
            for ligand in self._ligands_file_management.get_target_ligands_list(experiment_id, target_id):
                ligand_id = LigandId(ligand.ligand_id)
                for result in self._results_file_management.get_jobs_list(experiment_id, target_id, ligand_id):
                    job_metadata = self._results_file_management.get_job_metadata(experiment_id, target_id, ligand_id,
                                                                                  JobId(result.job_id))

                    is_available = False

                    if job_metadata.docking_method == "umol":
                        is_available = self._results_file_management.check_umol_result_data_available(experiment_id,
                                                                                                      target_id,
                                                                                                      ligand_id,
                                                                                                      JobId(
                                                                                                          result.job_id))
                    elif job_metadata.docking_method == "diffdock":
                        is_available = self._results_file_management.check_diffdock_result_data_available(experiment_id,
                                                                                                          target_id,
                                                                                                          ligand_id,
                                                                                                          JobId(
                                                                                                              result.job_id))

                    if is_available:
                        results_list.append(JobMetaData(job_id=result.job_id,
                                                     target_id=target_id.value,
                                                     ligand_id=ligand_id.value,
                                                     experiment_id=experiment_id.value,
                                                     folding_method=result.folding_method,
                                                     docking_method=result.docking_method))

        return GetAllResultsListResponse(results_list=results_list)


class GetAllJobsListFeature:
    def __init__(self, results_file_management: ResultsFileManagement,
                 target_file_management: TargetsFileManagement,
                 ligand_file_management: LigandsFileManagement):
        self._results_file_management = results_file_management
        self._targets_file_management = target_file_management
        self._ligands_file_management = ligand_file_management

    def handle(self, request: GetAllJobsListRequest) -> GetAllJobsListResponse:
        experiment_id = ExperimentId(request.experiment_id)

        jobs_list = []

        for target in self._targets_file_management.get_targets_list(experiment_id):
            target_id = TargetId(target.target_id)
            for ligand in self._ligands_file_management.get_target_ligands_list(experiment_id, target_id):
                ligand_id = LigandId(ligand.ligand_id)
                for job in self._results_file_management.get_jobs_list(experiment_id, target_id, ligand_id):
                    job_metadata = self._results_file_management.get_job_metadata(experiment_id, target_id, ligand_id,
                                                                                  JobId(job.job_id))

                    is_available = False

                    if job_metadata.docking_method == "umol":
                        is_available = self._results_file_management.check_umol_result_data_available(experiment_id,
                                                                                                      target_id,
                                                                                                      ligand_id,
                                                                                                      JobId(
                                                                                                          job.job_id))
                    elif job_metadata.docking_method == "diffdock":
                        is_available = self._results_file_management.check_diffdock_result_data_available(experiment_id,
                                                                                                          target_id,
                                                                                                          ligand_id,
                                                                                                          JobId(
                                                                                                              job.job_id))
                    if not is_available:
                        jobs_list.append(JobMetaData(job_id=job.job_id,
                                                     target_id=target_id.value,
                                                     ligand_id=ligand_id.value,
                                                     experiment_id=experiment_id.value,
                                                     folding_method=job.folding_method,
                                                     docking_method=job.docking_method))

        return GetAllJobsListResponse(jobs_list=jobs_list)


class GetResultsListForTargetLigandFeature:
    def __init__(self, results_file_management: ResultsFileManagement):
        self._results_file_management = results_file_management

    def handle(self, request: GetResultsListForTargetLigandRequest) -> GetResultsListForTargetLigandResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        results_list = []

        for result in self._results_file_management.get_jobs_list(experiment_id, target_id, ligand_id):
            job_metadata = self._results_file_management.get_job_metadata(experiment_id, target_id, ligand_id,
                                                                          JobId(result.job_id))

            is_available = False

            if job_metadata.docking_method == "umol":
                is_available = self._results_file_management.check_umol_result_data_available(experiment_id,
                                                                                              target_id,
                                                                                              ligand_id,
                                                                                              JobId(result.job_id))
            elif job_metadata.docking_method == "diffdock":
                is_available = self._results_file_management.check_diffdock_result_data_available(experiment_id,
                                                                                                  target_id,
                                                                                                  ligand_id,
                                                                                                  JobId(result.job_id))

            if is_available:
                results_list.append(JobMetaData(job_id=result.job_id,
                                                     target_id=target_id.value,
                                                     ligand_id=ligand_id.value,
                                                     experiment_id=experiment_id.value,
                                                     folding_method=result.folding_method,
                                                     docking_method=result.docking_method))

        return GetResultsListForTargetLigandResponse(results_list=results_list)


class GetJobsListForTargetLigandFeature:
    def __init__(self, results_file_management: ResultsFileManagement):
        self._results_file_management = results_file_management

    def handle(self, request: GetJobsListForTargetLigandRequest) -> GetJobsListForTargetLigandResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)

        jobs_list = []

        for job in self._results_file_management.get_jobs_list(experiment_id, target_id, ligand_id):

            job_metadata = self._results_file_management.get_job_metadata(experiment_id, target_id, ligand_id, JobId(job.job_id))

            is_available = False

            if job_metadata.docking_method == "umol":
                is_available = self._results_file_management.check_umol_result_data_available(experiment_id,
                                                                                         target_id,
                                                                                         ligand_id,
                                                                                         JobId(job.job_id))
            elif job_metadata.docking_method == "diffdock":
                is_available = self._results_file_management.check_diffdock_result_data_available(experiment_id,
                                                                                              target_id,
                                                                                              ligand_id,
                                                                                              JobId(job.job_id))

            if not is_available:
                jobs_list.append(JobMetaData(job_id=job.job_id,
                                             target_id=target_id.value,
                                             ligand_id=ligand_id.value,
                                             experiment_id=experiment_id.value,
                                             folding_method=job.folding_method,
                                             docking_method=job.docking_method))

        return GetJobsListForTargetLigandResponse(jobs_list=jobs_list)


class CheckResultDataAvailableFeature:
    def __init__(self, file_management: ResultsFileManagement):
        self._file_management = file_management

    def handle(self, request: CheckResultDataAvailableRequest) -> CheckResultDataAvailableResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        job_metadata = self._file_management.get_job_metadata(experiment_id, target_id, ligand_id,
                                                                      job_id)

        is_available = False
        if job_metadata.docking_method == "umol":
            is_available = self._file_management.check_umol_result_data_available(experiment_id,
                                                                                          target_id,
                                                                                          ligand_id,
                                                                                          job_id)
        elif job_metadata.docking_method == "diffdock":
            is_available = self._file_management.check_diffdock_result_data_available(experiment_id,
                                                                                              target_id,
                                                                                              ligand_id,
                                                                                              job_id)

        return CheckResultDataAvailableResponse(result_available=is_available)


class CheckMsaDataAvailableFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: CheckMsaDataAvailableRequest) -> CheckMsaDataAvailableResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)

        is_available = self._file_management.check_msa_exists(experiment_id, target_id)

        return CheckMsaDataAvailableResponse(is_available=is_available)


class CheckBindingPocketDataAvailableFeature:
    def __init__(self, file_management: ResultsFileManagement):
        self._file_management = file_management

    def handle(self, request: CheckPocketDataAvailableRequest) -> CheckPocketDataAvailableResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        is_available = self._file_management.check_binding_pocket_exist(experiment_id, target_id, ligand_id, job_id)

        return CheckPocketDataAvailableResponse(is_available=is_available)


class GetBJobBindingPocketDataFeature:
    def __init__(self, file_management: ResultsFileManagement):
        self._file_management = file_management

    def handle(self, request: GetJobBindingPocketDataRequest) -> GetJobBindingPocketDataResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        ligand_id = LigandId(request.ligand_id)
        job_id = JobId(request.job_id)

        pocket_ids = self._file_management.get_result_input_pocketIds(experiment_id, target_id, ligand_id, job_id)

        return GetJobBindingPocketDataResponse(pocket_ids=pocket_ids)


class CheckFoldingDataAvailableFeature:
    def __init__(self, file_management: TargetsFileManagement):
        self._file_management = file_management

    def handle(self, request: CheckFoldingDataAvailableRequest) -> CheckFoldingDataAvailableResponse:
        experiment_id = ExperimentId(request.experiment_id)
        target_id = TargetId(request.target_id)
        folding_method = request.folding_method

        is_available = self._file_management.check_folding_exist(experiment_id, target_id, folding_method)

        return CheckFoldingDataAvailableResponse(is_available=is_available)
