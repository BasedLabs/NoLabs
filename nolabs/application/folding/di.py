from typing import Annotated

from application.folding.use_cases import (GetJobFeature, GetJobStatusFeature,
                                           RunJobFeature, SetupJobFeature)
from esmfold_microservice import DefaultApi
from fastapi import Depends

from nolabs.infrastructure.di import InfrastructureDependencies


class FoldingDependencies:
    @staticmethod
    def run_job(
        esmfold: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.esmfold_microservice)
        ],
        esmfold_light: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.esmfold_light_microservice)
        ],
        rosettafold: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.esmfold_light_microservice)
        ],
    ) -> RunJobFeature:
        return RunJobFeature(
            esmfold=esmfold, esmfold_light=esmfold_light, rosettafold=rosettafold
        )

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()

    @staticmethod
    def get_job_status(
        esmfold: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.esmfold_microservice)
        ],
        esmfold_light: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.esmfold_light_microservice)
        ],
        rosettafold: Annotated[
            DefaultApi, Depends(InfrastructureDependencies.esmfold_light_microservice)
        ],
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            esmfold=esmfold, esmfold_light=esmfold_light, rosettafold=rosettafold
        )

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()
