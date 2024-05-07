__all__ = [
    'ConformationsDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.refined.application.use_cases.conformations.use_cases import SetupJobFeature, GetJobFeature, RunJobFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies


class ConformationsDependencies:
    @staticmethod
    def run_job(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.conformations_microservice)]
    ) -> RunJobFeature:
        return RunJobFeature(api=api)

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
