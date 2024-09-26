from nolabs.application.small_molecules_design.use_cases import (
    DeleteJobFeature,
    GetJobFeature,
    GetJobLogsFeature,
    GetJobSmilesFeature,
    GetJobStatusFeature,
    RunLearningStageJobFeature,
    RunSamplingStageJobFeature,
    SetupJobFeature,
    StopJobFeature,
)


class SmallMoleculesDesignDependencies:
    @staticmethod
    def delete_job() -> DeleteJobFeature:
        return DeleteJobFeature()

    @staticmethod
    def get_job_status() -> GetJobStatusFeature:
        return GetJobStatusFeature()

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def get_job_logs() -> GetJobLogsFeature:
        return GetJobLogsFeature()

    @staticmethod
    def get_job_smiles() -> GetJobSmilesFeature:
        return GetJobSmilesFeature()

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()

    @staticmethod
    def run_learning() -> RunLearningStageJobFeature:
        return RunLearningStageJobFeature()

    @staticmethod
    def run_sampling() -> RunSamplingStageJobFeature:
        return RunSamplingStageJobFeature()

    @staticmethod
    def stop_job() -> StopJobFeature:
        return StopJobFeature()
