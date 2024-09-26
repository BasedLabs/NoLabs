from nolabs.application.diffdock.use_cases import (
    GetJobFeature,
    RunJobFeature,
    SetupJobFeature,
)


class DiffDockDependencies:
    @staticmethod
    def run_job() -> RunJobFeature:
        return RunJobFeature()

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
