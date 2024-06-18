__all__ = [
    'GeneOntologyDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.application.use_cases.gene_ontology.use_cases import SetupJobFeature, RunJobFeature, GetJobFeature
from nolabs.infrastructure.di import InfrastructureDependencies


class GeneOntologyDependencies:
    @staticmethod
    def run_job(
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.gene_ontology_microservice)]
    ) -> RunJobFeature:
        return RunJobFeature(api=api)

    @staticmethod
    def setup_job() -> SetupJobFeature:
        return SetupJobFeature()

    @staticmethod
    def get_job() -> GetJobFeature:
        return GetJobFeature()