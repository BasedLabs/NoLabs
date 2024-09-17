__all__ = ["MsaGenerationDependencies"]

from typing import Annotated

import msa_light_microservice
from fastapi import Depends

from nolabs.application.use_cases.msa_generation.use_cases import (
    GetJobFeature, GetJobStatusFeature, RunJobFeature, SetupJobFeature)
from nolabs.infrastructure.di import InfrastructureDependencies


class MsaGenerationDependencies:
    @staticmethod
    def run_job(
        api: Annotated[
            msa_light_microservice.DefaultApi,
            Depends(InfrastructureDependencies.msa_light_microservice),
        ],
    ) -> RunJobFeature:
        return RunJobFeature(api=api)

    @staticmethod
    def get_job_status(
        api: Annotated[
            msa_light_microservice.DefaultApi,
            Depends(InfrastructureDependencies.msa_light_microservice),
        ]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(api=api)

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
