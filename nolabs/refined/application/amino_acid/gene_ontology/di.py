__all__ = [
    'GeneOntologyDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.refined.application.amino_acid.gene_ontology.use_cases import RunGeneOntologyFeature, GetGeneOntologyJobFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies
from nolabs.refined.infrastructure.settings import Settings


class GeneOntologyDependencies:
    @staticmethod
    def start(
            settings: Annotated[Settings, Depends(InfrastructureDependencies.settings)],
            api: Annotated[DefaultApi, Depends(InfrastructureDependencies.gene_ontology_microservice)]
    ) -> RunGeneOntologyFeature:
        return RunGeneOntologyFeature(settings=settings, api=api)

    @staticmethod
    def get_job() -> GetGeneOntologyJobFeature:
        return GetGeneOntologyJobFeature()
