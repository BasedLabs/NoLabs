__all__ = ["ProteinDesignDependencies"]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.application.use_cases.protein_design.use_cases import (
    GetJobFeature, GetJobStatusFeature, RunJobFeature, SetupJobFeature)
from nolabs.infrastructure.di import InfrastructureDependencies


class ProteinDesignDependencies:
    @staticmethod
    def run_job(
        api: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.protein_design_microservice)
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
            DefaultApi, Depends(InfrastructureDependencies.protein_design_microservice)
        ]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(api=api)
