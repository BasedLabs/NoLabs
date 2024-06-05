__all__ = [
    'CheckBioBuddyEnabledFeature'
]

import os
from uuid import UUID

import p2rank_microservice
from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.binding_pockets.api_models import GetJobStatusResponse, JobResponse, \
    SetupJobRequest
from nolabs.refined.application.use_cases.biobuddy.api_models import CheckBioBuddyEnabledResponse
from nolabs.refined.domain.models.common import Protein, JobId, JobName, Experiment
from nolabs.refined.domain.models.pocket_prediction import PocketPredictionJob
from nolabs.utils import generate_uuid

class CheckBioBuddyEnabledFeature:
    def __init__(self):
        pass

    def handle(self) -> CheckBioBuddyEnabledResponse:
        enable_biobuddy = os.getenv('ENABLE_BIOBUDDY', 'false').lower() == 'true'
        return CheckBioBuddyEnabledResponse(enabled=enable_biobuddy)

