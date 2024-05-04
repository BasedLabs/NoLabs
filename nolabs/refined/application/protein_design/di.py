__all__ = [
    'ProteinDesignDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.refined.application.amino_acid.localisation.use_cases import RunLocalisationFeature, GetLocalisationJobFeature
from nolabs.refined.application.protein_design.use_cases import RunProteinDesignFeature, GetProteinDesignJobFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.refined.infrastructure.settings import Settings


class ProteinDesignDependencies:
    @staticmethod
    def start(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.protein_design_microservice)]
    ) -> RunProteinDesignFeature:
        return RunProteinDesignFeature(api=api)

    @staticmethod
    def get_job() -> GetProteinDesignJobFeature:
        return GetProteinDesignJobFeature()
