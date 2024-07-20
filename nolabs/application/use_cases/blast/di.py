from typing import Annotated

from blast_query_microservice import DefaultApi
from fastapi import Depends

from nolabs.application.use_cases.blast.use_cases import GetJobFeature, RunJobFeature, \
    GetJobStatusFeature, SetupJobFeature
from nolabs.infrastructure.di import InfrastructureDependencies


class BlastDependencies:
    @staticmethod
    def run_job(
            blast: Annotated[DefaultApi, Depends(InfrastructureDependencies.blast_query_microservice)]
    ) -> RunJobFeature:
        return RunJobFeature(
            blast=blast
        )

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def get_job_status(blast: Annotated[
        DefaultApi, Depends(InfrastructureDependencies.blast_query_microservice)]) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            blast=blast
        )

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
