from typing import Annotated

from fastapi import APIRouter, Depends

from nolabs.refined.application.use_cases.biobuddy.api_models import CheckBioBuddyEnabledResponse
from nolabs.refined.application.use_cases.biobuddy.di import BiobuddyDependencies
from nolabs.refined.application.use_cases.biobuddy.use_cases import CheckBioBuddyEnabledFeature

router = APIRouter(
    prefix='/api/v1/biobuddy',
    tags=['biobuddy']
)

@router.get('/check_biobuddy_enabled')
async def check_biobuddy_enabled(
        feature: Annotated[
    CheckBioBuddyEnabledFeature, Depends(BiobuddyDependencies.check_biobuddy_enabled)]) -> CheckBioBuddyEnabledResponse:
    return feature.handle()