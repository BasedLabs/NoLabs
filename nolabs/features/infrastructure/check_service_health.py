import urllib.parse

import aiohttp

from nolabs.infrastructure.settings import Settings

from nolabs.api_models.drug_discovery import CheckServiceHealthyResponse


async def _healthcheck(host: str, path: str = 'livez') -> CheckServiceHealthyResponse:
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(urllib.parse.urljoin(host, path)) as resp:
                return CheckServiceHealthyResponse(is_healthy=(resp.status // 100 == 2))
    except aiohttp.ClientError as e:
        return CheckServiceHealthyResponse(is_healthy=False)


class CheckUmolServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self) -> CheckServiceHealthyResponse:
        return await _healthcheck(self._settings.umol_host, 'docs')


class CheckDiffDockServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self) -> CheckServiceHealthyResponse:
        return await _healthcheck(self._settings.diffdock_host, 'docs')


class CheckEsmFoldServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self) -> CheckServiceHealthyResponse:
        return await _healthcheck(self._settings.esmfold_host, 'docs')


class CheckEsmFoldLightServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self) -> CheckServiceHealthyResponse:
        return await _healthcheck(self._settings.esmfold_light_host, 'docs')


class CheckMsaServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self) -> CheckServiceHealthyResponse:
        return await _healthcheck(self._settings.msa_light_host, 'docs')


class CheckP2RankServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self) -> CheckServiceHealthyResponse:
        return await _healthcheck(self._settings.p2rank_host, 'docs')


class CheckRosettaFoldServiceHealthFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    async def handle(self) -> CheckServiceHealthyResponse:
        return await _healthcheck(self._settings.rosettafold_host)
