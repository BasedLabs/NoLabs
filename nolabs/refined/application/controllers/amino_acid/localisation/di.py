__all__ = [
    'LocalisationDependencies'
]


from typing import Annotated

from fastapi import Depends

from nolabs.refined.application.usecases.amino_acid.localisation.use_cases import RunLocalisationFeature, \
    GetJobsMetadataFeature, GetJobFeature, DeleteJobFeature, UpdateJobFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.refined.infrastructure.settings import Settings


class LocalisationDependencies:
    @staticmethod
    def start(
            settings: Annotated[Settings, Depends(InfrastructureDependencies.settings)]
    ) -> RunLocalisationFeature:
        return RunLocalisationFeature(settings=settings)

    @staticmethod
    def job_metadata() -> GetJobsMetadataFeature:
        return GetJobsMetadataFeature()

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def delete_job():
        return DeleteJobFeature()

    @staticmethod
    def update_job():
        return UpdateJobFeature()
