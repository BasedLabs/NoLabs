from nolabs.features.drug_discovery.data_models.result import JobId
from nolabs.infrastructure.settings import Settings

# TODO: Add self-hosted MSA microservice
# if settings.drug_discovery_self_hosted_msa:
from msa_light_microservice import DefaultApi as MsaDefaultApi
from msa_light_microservice import ApiClient as MsaApiClient
from msa_light_microservice import Configuration as MsaConfiguration

from p2rank_microservice import DefaultApi as P2RankDefaultApi
from p2rank_microservice import ApiClient as P2RankApiClient
from p2rank_microservice import Configuration as P2RankConfiguration

from umol_microservice import DefaultApi as UmolDefaultApi
from umol_microservice import ApiClient as UmolApiClient
from umol_microservice import Configuration as UmolConfiguration

settings = Settings()
if settings.is_light_infrastructure:
    from esmfold_light_microservice import DefaultApi as FoldingDefaultApi
    from esmfold_light_microservice import ApiClient as FoldingApiClient
    from esmfold_light_microservice import Configuration as FoldingConfiguration
else:
    from esmfold_microservice import DefaultApi as FoldingDefaultApi
    from esmfold_microservice import ApiClient as FoldingApiClient
    from esmfold_microservice import Configuration as FoldingConfiguration


from nolabs.api_models.drug_discovery import CheckJobIsRunningRequest, CheckJobIsRunningResponse

class CheckUmolRunningFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self, request: CheckJobIsRunningRequest) -> CheckJobIsRunningResponse:
        assert request

        job_id = JobId(request.job_id)

        configuration = UmolConfiguration(
            host=self._settings.umol_host,
        )
        with UmolApiClient(configuration=configuration) as client:
            api_instance = UmolDefaultApi(client)
            response = api_instance.is_job_running_job_job_id_is_running_get(job_id=job_id.value).is_running

        return CheckJobIsRunningResponse(is_running=response)

class CheckFoldingRunningFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self, request: CheckJobIsRunningRequest) -> CheckJobIsRunningResponse:
        assert request

        job_id = JobId(request.job_id)

        configuration = FoldingConfiguration(
            host=self._settings.umol_host,
        )
        with FoldingApiClient(configuration=configuration) as client:
            api_instance = FoldingDefaultApi(client)
            response = api_instance.is_job_running_job_job_id_is_running_get(job_id=job_id.value).is_running

        return CheckJobIsRunningResponse(is_running=response)



class CheckMsaRunningFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self, request: CheckJobIsRunningRequest) -> CheckJobIsRunningResponse:
        job_id = JobId(request.job_id)

        configuration = MsaConfiguration(
            host=self._settings.msa_light_host,
        )
        with MsaApiClient(configuration=configuration) as client:
            api_instance = MsaDefaultApi(client)
            response = api_instance.is_job_running_job_job_id_is_running_get(job_id=job_id.value).is_running

        return CheckJobIsRunningResponse(is_running=response)


class CheckP2RankRunningFeature:
    def __init__(self, settings: Settings):
        self._settings = settings

    def handle(self, request: CheckJobIsRunningRequest) -> CheckJobIsRunningResponse:
        job_id = JobId(request.job_id)

        configuration = P2RankConfiguration(
            host=self._settings.p2rank_host,
        )
        with P2RankApiClient(configuration=configuration) as client:
            api_instance = P2RankDefaultApi(client)
            response = api_instance.is_job_running_job_job_id_is_running_get(job_id=job_id.value).is_running

        return CheckJobIsRunningResponse(is_running=response)
