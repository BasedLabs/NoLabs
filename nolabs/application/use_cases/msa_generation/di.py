__all__ = ['MsaGenerationDependencies']

from typing import Annotated

import msa_light_microservice
from fastapi import Depends

from nolabs.application.use_cases.msa_generation.use_cases import RunJobFeature, GetJobStatusFeature, GetJobFeature, \
    SetupJobFeature
from nolabs.infrastructure.settings import MsaLightMicroserviceSettings
from nolabs.infrastructure.di import InfrastructureDependencies


class MsaGenerationDependencies:
    @staticmethod
    def run_job(
            settings: Annotated[MsaLightMicroserviceSettings, Depends(InfrastructureDependencies.msa_light_settings)],
            api: Annotated[msa_light_microservice.DefaultApi, Depends(InfrastructureDependencies.msa_light_microservice)],
    ) -> RunJobFeature:
        return RunJobFeature(
            settings=settings,
            api=api
        )

    @staticmethod
    def get_job_status(
            api: Annotated[msa_light_microservice.DefaultApi, Depends(InfrastructureDependencies.msa_light_microservice)]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            api=api
        )

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()