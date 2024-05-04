__all__ = [
    'ConformationsDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.refined.application.conformations.use_cases import GetConformationsJobFeature, RunConformationsJobFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies


class ConformationsDependencies:
    @staticmethod
    def start(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.conformations_microservice)]
    ) -> RunConformationsJobFeature:
        return RunConformationsJobFeature(api=api)

    @staticmethod
    def get_job() -> GetConformationsJobFeature:
        return GetConformationsJobFeature()
