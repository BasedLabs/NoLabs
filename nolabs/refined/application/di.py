__all__ = [
    'EventHandlersContainer',
    'ApplicationContainer'
]


from dependency_injector import providers, containers

import nolabs.refined.application.features.experiments as experiments
import nolabs.refined.application.event_handlers as event_handlers
import nolabs.refined.application.features.amino_acid.localisation as localisation


# region Features

class ExperimentsFeaturesContainer(containers.DeclarativeContainer):
    get_experiments_metadata_feature = providers.Factory(experiments.GetExperimentsMetadataFeature)
    delete_experiment_feature = providers.Factory(experiments.DeleteExperimentFeature)
    change_experiment_name_feature = providers.Factory(experiments.ChangeExperimentNameFeature)
    create_experiment_feature = providers.Factory(experiments.CreateExperimentFeature)


class AminoAcidFeaturesContainer(containers.DeclarativeContainer):
    # region Localisation

    get_jobs_feature = providers.Factory(localisation.GetJobsMetadataFeature)
    get_job_feature = providers.Factory(localisation.GetJobFeature)
    run_localisation_feature = providers.Factory(localisation.RunLocalisationFeature)
    delete_job_feature = providers.Factory(localisation.DeleteJobFeature)
    update = providers.Factory(localisation.UpdateJobFeature)

    # endregion

# endregion


class EventHandlersContainer(containers.DeclarativeContainer):
    __self__ = providers.Self()

    aa_created = providers.Factory(event_handlers.AminoAcidCreatedEventHandler)
    protein_created = providers.Factory(event_handlers.ProteinCreatedEventHandler)


class ApplicationContainer(containers.DeclarativeContainer):
    experiments_features = providers.Container(
        ExperimentsFeaturesContainer
    )
    amino_acid_features = providers.Container(
        AminoAcidFeaturesContainer
    )
    event_handlers = providers.Container(EventHandlersContainer)
