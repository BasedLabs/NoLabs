from typing import Annotated

from esmfold_microservice import DefaultApi
from fastapi import Depends

from nolabs.refined.application.amino_acid.folding.use_cases import GetFoldingJobFeature, RunFoldingFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.refined.infrastructure.settings import Settings


class FoldingDependencies:
    @staticmethod
    def start(
            settings: Annotated[Settings, Depends(InfrastructureDependencies.settings)],
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.localisation_microservice)]
    ) -> RunFoldingFeature:
        return RunFoldingFeature(settings=settings, api=api)

    @staticmethod
    def get_job() -> GetFoldingJobFeature:
        return GetFoldingJobFeature()