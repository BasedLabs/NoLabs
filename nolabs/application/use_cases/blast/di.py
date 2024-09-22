from nolabs.application.use_cases.blast.use_cases import (GetJobFeature,
                                                          RunJobFeature,
                                                          SetupJobFeature)


class BlastDependencies:
    @staticmethod
    def run_job() -> RunJobFeature:
        return RunJobFeature()

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
