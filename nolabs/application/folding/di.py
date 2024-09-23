from nolabs.application.folding.use_cases import RunJobFeature, GetJobFeature, SetupJobFeature


class FoldingDependencies:
    @staticmethod
    def run_job() -> RunJobFeature:
        return RunJobFeature()

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
