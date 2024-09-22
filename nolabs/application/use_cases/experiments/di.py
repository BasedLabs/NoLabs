__all__ = ["ExperimentsDependencies"]


from nolabs.application.use_cases.experiments.use_cases import (
    CreateExperimentFeature, DeleteExperimentFeature,
    GetExperimentsMetadataFeature, UpdateExperimentFeature)


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
