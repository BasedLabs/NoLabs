from requests import RequestException

from nolabs.infrastructure.settings import Settings
import requests

from nolabs.api_models.drug_discovery import CheckServiceHealthyResponse


class CheckUmolServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self) -> CheckServiceHealthyResponse:

        host = self._settings.umol_host
        # Attempt to ping the host briefly to check if it's healthy
        try:
            response = requests.get(f"{host}/docs", timeout=0.5)  # Timeout set to 5 seconds
            # A 2xx status code indicates success
            is_healthy = response.status_code // 100 == 2
        except RequestException:
            # Any exception (connection error, timeout, etc.) implies the service is not healthy
            is_healthy = False

        return CheckServiceHealthyResponse(is_healthy=is_healthy)


class CheckFoldingServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self) -> CheckServiceHealthyResponse:

        host = None
        if self._settings.is_light_infrastructure:
            host = self._settings.esmfold_light_host
        else:
            host = self._settings.esmfold_host

        try:
            response = requests.get(f"{host}/docs", timeout=0.5)  # Timeout set to 0.1 seconds
            # A 2xx status code indicates success
            is_healthy = response.status_code // 100 == 2
        except RequestException:
            # Any exception (connection error, timeout, etc.) implies the service is not healthy
            is_healthy = False

        return CheckServiceHealthyResponse(is_healthy=is_healthy)


class CheckMsaServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self) -> CheckServiceHealthyResponse:

        host = self._settings.msa_light_host

        try:
            response = requests.get(f"{host}/docs", timeout=0.5)  # Timeout set to 0.1 seconds
            # A 2xx status code indicates success
            is_healthy = response.status_code // 100 == 2
        except RequestException:
            # Any exception (connection error, timeout, etc.) implies the service is not healthy
            is_healthy = False

        return CheckServiceHealthyResponse(is_healthy=is_healthy)


class CheckP2RankServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self) -> CheckServiceHealthyResponse:
        host = self._settings.p2rank_host

        try:
            response = requests.get(f"{host}/docs", timeout=0.5)  # Timeout set to 0.1 seconds
            # A 2xx status code indicates success
            is_healthy = response.status_code // 100 == 2
        except RequestException:
            # Any exception (connection error, timeout, etc.) implies the service is not healthy
            is_healthy = False

        return CheckServiceHealthyResponse(is_healthy=is_healthy)
