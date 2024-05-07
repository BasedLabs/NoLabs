__all__ = [
    'ExperimentsDependencies'
]


from nolabs.refined.application.use_cases.experiments.use_cases import GetExperimentsMetadataFeature, \
    CreateExperimentFeature, DeleteExperimentFeature, UpdateExperimentFeature


class ExperimentsDependencies:
    @staticmethod
    def experiments() -> GetExperimentsMetadataFeature:
        return GetExperimentsMetadataFeature()

    @staticmethod
    def create_experiment() -> CreateExperimentFeature:
        return CreateExperimentFeature()

    @staticmethod
    def delete_experiment() -> DeleteExperimentFeature:
        return DeleteExperimentFeature()

    @staticmethod
    def update_experiment() -> UpdateExperimentFeature:
        return UpdateExperimentFeature()
