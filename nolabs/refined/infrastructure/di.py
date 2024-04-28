from dependency_injector import containers, providers

from nolabs.refined.domain.repository import Repository
from nolabs.refined.infrastructure.settings import Settings


class InfrastructureContainer(containers.DeclarativeContainer):
    settings = providers.Singleton(Settings)
    repository = providers.Singleton(Repository)
