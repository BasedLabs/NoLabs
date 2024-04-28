from dependency_injector import containers, providers

from nolabs.refined.application.features.experiments import GetExperimentsMetadataFeature, DeleteExperimentFeature, \
    ChangeExperimentNameFeature, CreateExperimentFeature


class ExperimentsContainer(containers.DeclarativeContainer):
