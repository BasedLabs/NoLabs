import os

from nolabs.api_models.biobuddy import CheckBioBuddyEnabledResponse


class CheckBioBuddyEnabledFeature:
    def __init__(self):
        pass

    def handle(self) -> CheckBioBuddyEnabledResponse:
        enable_biobuddy = os.getenv('ENABLE_BIOBUDDY', 'false').lower() == 'true'
        return CheckBioBuddyEnabledResponse(enabled=enable_biobuddy)

