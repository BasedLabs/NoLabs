__all__ = ["LocalisationDependencies"]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.application.use_cases.localisation.use_cases import (
    GetJobFeature, GetJobStatusFeature, RunJobFeature, SetupJobFeature)
from nolabs.infrastructure.di import InfrastructureDependencies


class LocalisationDependencies:
    @staticmethod
    def run_job(
        api: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.localisation_microservice)
        ]
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
            DefaultApi, Depends(InfrastructureDependencies.localisation_microservice)
        ]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(api=api)
