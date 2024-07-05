__all__ = [
    'BindingPocketsDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.application.use_cases.binding_pockets.use_cases import RunJobFeature, GetJobFeature, GetJobStatusFeature, \
    SetupJobFeature
from nolabs.infrastructure.di import InfrastructureDependencies


class BindingPocketsDependencies:
    @staticmethod
    def run_job(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.p2rank_microservice)]
    ) -> RunJobFeature:
        return RunJobFeature(api=api)

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def get_job_status(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.p2rank_microservice)]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            api=api
        )

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
