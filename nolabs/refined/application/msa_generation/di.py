__all__ = ['MsaGenerationDependencies']

from typing import Annotated

import msa_light_microservice
from fastapi import Depends

from nolabs.refined.application.msa_generation.use_cases import RunMsaGenerationJobFeature, GetJobStatusFeature, \
    GetMsaPredictionJobResultFeature
from nolabs.refined.infrastructure.settings import MsaLightMicroserviceSettings
from nolabs.refined.infrastructure.di import InfrastructureDependencies


class MsaGenerationDependencies:
    @staticmethod
    def start(
            settings: Annotated[MsaLightMicroserviceSettings, Depends(InfrastructureDependencies.msa_light_settings)],
            api: Annotated[msa_light_microservice.DefaultApi, Depends(InfrastructureDependencies.msa_light_microservice)],
    ) -> RunMsaGenerationJobFeature:
        return RunMsaGenerationJobFeature(
            settings=settings,
            api=api
        )

    @staticmethod
    def get_status(
            api: Annotated[msa_light_microservice.DefaultApi, Depends(InfrastructureDependencies.msa_light_microservice)]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            api=api
        )

    @staticmethod
    def get_job() -> GetMsaPredictionJobResultFeature:
        return GetMsaPredictionJobResultFeature()