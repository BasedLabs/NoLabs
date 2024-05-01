__all__ = [
    'SolubilityDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.refined.application.amino_acid.localisation.use_cases import RunLocalisationFeature, GetLocalisationJobFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.refined.infrastructure.settings import Settings


class SolubilityDependencies:
    @staticmethod
    def start(
            settings: Annotated[Settings, Depends(InfrastructureDependencies.settings)],
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.localisation_microservice)]
    ) -> RunLocalisationFeature:
        return RunLocalisationFeature(settings=settings, api=api)

    @staticmethod
    def get_job() -> GetLocalisationJobFeature:
        return GetLocalisationJobFeature()
