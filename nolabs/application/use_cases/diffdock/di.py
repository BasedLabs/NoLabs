from typing import Annotated

from diffdock_microservice import DefaultApi
from fastapi import Depends

from nolabs.application.use_cases.diffdock.use_cases import GetJobFeature, RunJobFeature, \
    GetJobStatusFeature, SetupJobFeature
from nolabs.infrastructure.di import InfrastructureDependencies


class DiffDockDependencies:
    @staticmethod
    def run_job(
            diffdock: Annotated[DefaultApi, Depends(InfrastructureDependencies.diffdock_microservice)]
    ) -> RunJobFeature:
        return RunJobFeature(
            diffdock=diffdock
        )

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def get_job_status(diffdock: Annotated[
        DefaultApi, Depends(InfrastructureDependencies.diffdock_microservice)]) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            diffdock=diffdock
        )

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
