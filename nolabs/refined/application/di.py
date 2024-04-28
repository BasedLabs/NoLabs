from dependency_injector import containers, providers

from nolabs.refined.application.event_handlers.amino_acid_event_handlers import AminoAcidCreatedEventHandler
from nolabs.refined.application.event_handlers.protein_event_handlers import ProteinCreatedEventHandler
from nolabs.refined.application.features.experiments import GetExperimentsMetadataFeature, DeleteExperimentFeature, \
    ChangeExperimentNameFeature, CreateExperimentFeature
from nolabs.refined.settings import Settings


class ExperimentsControllerContainer(containers.DeclarativeContainer):
    get_experiments_metadata_feature = providers.Factory(GetExperimentsMetadataFeature)
    delete_experiment_feature = providers.Factory(DeleteExperimentFeature)
    change_experiment_name_feature = providers.Factory(ChangeExperimentNameFeature)
    create_experiment_feature = providers.Factory(CreateExperimentFeature)


class EventHandlersContainer(containers.DeclarativeContainer):
    __self__ = providers.Self()

    aa_created = providers.Factory(AminoAcidCreatedEventHandler)
    protein_created = providers.Factory(ProteinCreatedEventHandler)


class Container(containers.DeclarativeContainer):
    settings = providers.Singleton(Settings)

    # event handlers
