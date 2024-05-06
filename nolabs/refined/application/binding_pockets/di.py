__all__ = [
    'BindingPocketsDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.refined.application.binding_pockets.use_cases import GetBindingPocketsFeature, PredictBindingPocketsFeature, \
    SetManualBindingPocketsFeature, GetJobStatusFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies


class BindingPocketsDependencies:
    @staticmethod
    def predict_binding_pockets(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.p2rank_microservice)]
    ) -> PredictBindingPocketsFeature:
        return PredictBindingPocketsFeature(api=api)

    @staticmethod
    def get_binding_pockets() -> GetBindingPocketsFeature:
        return GetBindingPocketsFeature()

    @staticmethod
    def get_status(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.p2rank_microservice)]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(
            api=api
        )

    @staticmethod
    def set_binding_pockets() -> SetManualBindingPocketsFeature:
        return SetManualBindingPocketsFeature()
