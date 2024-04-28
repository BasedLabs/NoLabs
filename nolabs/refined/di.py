from dependency_injector import containers, providers

import nolabs.refined.application.di as di_application
import nolabs.refined.infrastructure.di as di_infrastructure


class Nolabs(containers.DeclarativeContainer):
    application = providers.Container(
        di_application.ApplicationContainer
    )

    infrastructure = providers.Container(
        di_infrastructure.InfrastructureContainer
    )
