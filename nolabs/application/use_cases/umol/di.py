from typing import Annotated

from diffdock_microservice import DefaultApi
from fastapi import Depends

from nolabs.application.use_cases.umol.use_cases import GetJobFeature, RunJobFeature, \
    GetJobStatusFeature, SetupJobFeature
from nolabs.infrastructure.di import InfrastructureDependencies


class UmolDependencies:
    @staticmethod
    def run_job(
            umol: Annotated[DefaultApi, Depends(InfrastructureDependencies.umol_microservice)]
    ) -> RunJobFeature:
        return RunJobFeature(
            umol=umol
        )

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def get_job_status(
            umol: Annotated[DefaultApi, Depends(InfrastructureDependencies.umol_microservice)]) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            umol=umol
        )

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
