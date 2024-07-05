__all__ = [
    'SolubilityDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.application.use_cases.solubility.use_cases import RunJobFeature, GetJobFeature, SetupJobFeature, \
    GetJobStatusFeature
from nolabs.infrastructure.di import InfrastructureDependencies


class SolubilityDependencies:
    @staticmethod
    def run_job(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.solubility_microservice)]
    ) -> RunJobFeature:
        return RunJobFeature(api=api)

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()

    @staticmethod
    def get_job_status(
            api: Annotated[
                DefaultApi, Depends(InfrastructureDependencies.solubility_microservice)]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            api=api
        )