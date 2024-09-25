import biobuddy_microservice
import diffdock_microservice
import external_data_query_microservice
import reinvent_microservice

from nolabs.infrastructure.settings import settings


class InfrastructureDependencies:
    @staticmethod
    def biobuddy_microservice() -> biobuddy_microservice.DefaultApi:
        configuration = biobuddy_microservice.Configuration(host=settings.biobuddy_host)
        client = biobuddy_microservice.ApiClient(configuration=configuration)
        return biobuddy_microservice.DefaultApi(client)

    @staticmethod
    def external_query_microservice() -> external_data_query_microservice.DefaultApi:
        configuration = external_data_query_microservice.Configuration(
            host=settings.external_query_host
        )
        client = external_data_query_microservice.ApiClient(configuration=configuration)
        return external_data_query_microservice.DefaultApi(client)

    @staticmethod
    def reinvent_microservice() -> reinvent_microservice.ReinventApi:
        configuration = reinvent_microservice.Configuration(host=settings.reinvent_host)
        client = reinvent_microservice.ApiClient(configuration=configuration)
        return reinvent_microservice.ReinventApi(api_client=client)

    @staticmethod
    def diffdock_microservice() -> diffdock_microservice.DefaultApi:
        configuration = diffdock_microservice.Configuration(host=settings.diffdock_host)
        client = diffdock_microservice.ApiClient(configuration=configuration)
        return diffdock_microservice.DefaultApi(api_client=client)
