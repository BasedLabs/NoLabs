import biobuddy_microservice
import chembl_query_microservice
import localisation_microservice
import esmfold_microservice
import esmfold_light_microservice
import rcsb_pdb_query_microservice
import rosettafold_microservice
import gene_ontology_microservice
import reinvent_microservice
import solubility_microservice
import protein_design_microservice
import conformations_microservice
import p2rank_microservice
import msa_light_microservice
import diffdock_microservice
import umol_microservice

from nolabs.refined.infrastructure.settings import Settings, MsaLightMicroserviceSettings


class InfrastructureDependencies:
    @staticmethod
    def settings() -> Settings:
        return Settings.load()

    @staticmethod
    def biobuddy_microservice() -> biobuddy_microservice.DefaultApi:
        settings = Settings.load()
        configuration = biobuddy_microservice.Configuration(
            host=settings.biobuddy.microservice
        )
        client = biobuddy_microservice.ApiClient(configuration=configuration)
        return biobuddy_microservice.DefaultApi(client)

    @staticmethod
    def chembl_query_microservice() -> chembl_query_microservice.DefaultApi:
        settings = Settings.load()
        configuration = chembl_query_microservice.Configuration(
            host=settings.chembl.microservice
        )
        client = chembl_query_microservice.ApiClient(configuration=configuration)
        return chembl_query_microservice.DefaultApi(client)

    @staticmethod
    def rcsb_pdb_query_microservice() -> rcsb_pdb_query_microservice.DefaultApi:
        settings = Settings.load()
        configuration = rcsb_pdb_query_microservice.Configuration(
            host=settings.rcsb_pdb.microservice
        )
        client = rcsb_pdb_query_microservice.ApiClient(configuration=configuration)
        return rcsb_pdb_query_microservice.DefaultApi(client)


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
    def esmfold_light_microservice() -> esmfold_light_microservice.DefaultApi:
        settings = Settings.load()
        configuration = esmfold_light_microservice.Configuration(
            host=settings.esmfold_light.microservice
        )
        client = esmfold_light_microservice.ApiClient(configuration=configuration)
        return esmfold_light_microservice.DefaultApi(client)

    @staticmethod
    def rosettafold_microservice() -> rosettafold_microservice.DefaultApi:
        settings = Settings.load()
        configuration = rosettafold_microservice.Configuration(
            host=settings.rosettafold.microservice
        )
        client = rosettafold_microservice.ApiClient(configuration=configuration)
        return rosettafold_microservice.DefaultApi(client)

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
            host=settings.solubility.microservice
        )
        client = solubility_microservice.ApiClient(configuration=configuration)
        return solubility_microservice.DefaultApi(client)

    @staticmethod
    def reinvent_microservice() -> reinvent_microservice.ReinventApi:
        settings = Settings.load()
        configuration = reinvent_microservice.Configuration(
            host=settings.reinvent.microservice
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

    @staticmethod
    def p2rank_microservice() -> p2rank_microservice.DefaultApi:
        settings = Settings.load()
        configuration = p2rank_microservice.Configuration(
            host=settings.p2rank.microservice
        )
        client = p2rank_microservice.ApiClient(configuration=configuration)
        return p2rank_microservice.DefaultApi(client)

    @staticmethod
    def msa_light_microservice() -> msa_light_microservice.DefaultApi:
        settings = Settings.load()
        configuration = msa_light_microservice.Configuration(
            host=settings.msa_light_microservice.microservice
        )
        client = msa_light_microservice.ApiClient(configuration=configuration)
        return msa_light_microservice.DefaultApi(client)

    @staticmethod
    def msa_light_settings() -> MsaLightMicroserviceSettings:
        settings = Settings.load()
        return settings.msa_light

    @staticmethod
    def diffdock_microservice() -> diffdock_microservice.DefaultApi:
        settings = Settings.load()
        configuration = diffdock_microservice.Configuration(
            host=settings.diffdock.microservice
        )
        client = diffdock_microservice.ApiClient(configuration=configuration)
        return diffdock_microservice.DefaultApi(api_client=client)

    @staticmethod
    def umol_microservice() -> umol_microservice.DefaultApi:
        settings = Settings.load()
        configuration = umol_microservice.Configuration(
            host=settings.umol_microservice.microservice
        )
        client = umol_microservice.ApiClient(configuration=configuration)
        return umol_microservice.DefaultApi(api_client=client)