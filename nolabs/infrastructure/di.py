import biobuddy_microservice
import external_data_query_microservice
from biobuddy_microservice.api.default_api import DefaultApi as BBDefaultApi
from external_data_query_microservice.api.default_api import DefaultApi

from nolabs.infrastructure.settings import settings


class InfrastructureDependencies:
    @staticmethod
    def biobuddy_microservice() -> BBDefaultApi:
        configuration = biobuddy_microservice.Configuration(host=settings.biobuddy_host)
        client = biobuddy_microservice.ApiClient(configuration=configuration)
        return biobuddy_microservice.DefaultApi(client)

    @staticmethod
    def external_query_microservice() -> DefaultApi:
        configuration = external_data_query_microservice.Configuration(
            host=settings.external_query_host
        )
        client = external_data_query_microservice.ApiClient(configuration=configuration)
        return external_data_query_microservice.DefaultApi(client)
