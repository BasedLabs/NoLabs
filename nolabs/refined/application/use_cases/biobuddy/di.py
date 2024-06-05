__all__ = [
    'BiobuddyDependencies'
]

from typing import Annotated

from fastapi import Depends
from localisation_microservice import DefaultApi

from nolabs.refined.application.use_cases.binding_pockets.use_cases import RunJobFeature, GetJobFeature, GetJobStatusFeature, \
    SetupJobFeature
from nolabs.refined.application.use_cases.biobuddy.use_cases import CheckBioBuddyEnabledFeature
from nolabs.refined.infrastructure.di import InfrastructureDependencies


class BiobuddyDependencies:
    @staticmethod
    def check_biobuddy_enabled() -> CheckBioBuddyEnabledFeature:
        return CheckBioBuddyEnabledFeature()