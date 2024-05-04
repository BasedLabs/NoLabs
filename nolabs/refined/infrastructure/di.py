import localisation_microservice
import esmfold_microservice
import gene_ontology_microservice
import reinvent_microservice
import solubility_microservice
import protein_design_microservice
import conformations_microservice

from nolabs.refined.infrastructure.settings import Settings


class InfrastructureDependencies:
    @staticmethod
    def settings() -> Settings:
        return Settings.load()

    @staticmethod
    def localisation_microservice() -> localisation_microservice.DefaultApi:
        settings = Settings.load()
        configuration = localisation_microservice.Configuration(
            host=settings.localisation.microservice
        )
        client = localisation_microservice.ApiClient(configuration=configuration)
        return localisation_microservice.DefaultApi(client)

    @staticmethod
    def esmfold_microservice() -> esmfold_microservice.DefaultApi:
        settings = Settings.load()
        configuration = esmfold_microservice.Configuration(
            host=settings.esmfold.microservice
        )
        client = esmfold_microservice.ApiClient(configuration=configuration)
        return esmfold_microservice.DefaultApi(client)

    @staticmethod
    def gene_ontology_microservice() -> gene_ontology_microservice.DefaultApi:
        settings = Settings.load()
        configuration = gene_ontology_microservice.Configuration(
            host=settings.gene_ontology.microservice
        )
        client = gene_ontology_microservice.ApiClient(configuration=configuration)
        return gene_ontology_microservice.DefaultApi(client)

    @staticmethod
    def solubility_microservice() -> solubility_microservice.DefaultApi:
        settings = Settings.load()
        configuration = solubility_microservice.Configuration(
            host=settings.solubility_microservice.microservice
        )
        client = solubility_microservice.ApiClient(configuration=configuration)
        return solubility_microservice.DefaultApi(client)

    @staticmethod
    def reinvent_microservice() -> reinvent_microservice.ReinventApi:
        settings = Settings.load()
        configuration = reinvent_microservice.Configuration(
            host=settings.reinvent_microservice.microservice
        )
        client = reinvent_microservice.ApiClient(configuration=configuration)
        return reinvent_microservice.ReinventApi(api_client=client)

    @staticmethod
    def protein_design_microservice() -> protein_design_microservice.DefaultApi:
        settings = Settings.load()
        configuration = protein_design_microservice.Configuration(
            host=settings.protein_design.microservice
        )
        client = protein_design_microservice.ApiClient(configuration=configuration)
        return protein_design_microservice.DefaultApi(client)

    @staticmethod
    def conformations_microservice() -> conformations_microservice.DefaultApi:
        settings = Settings.load()
        configuration = protein_design_microservice.Configuration(
            host=settings.conformations.microservice
        )
        client = conformations_microservice.ApiClient(configuration=configuration)
        return conformations_microservice.DefaultApi(client)