__all__ = [
    'Settings'
]

import os

from pydantic import BaseModel
from pydantic_settings import BaseSettings

from nolabs.refined.infrastructure.environment import Environment


class MicroserviceSettings(BaseModel):
    microservice: str = ''


class MsaLightMicroserviceSettings(MicroserviceSettings):
    msa_server_url: str = ''


settings_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                             f'appsettings.{os.environ["NOLABS_ENVIRONMENT"]}.json')


class Settings(BaseSettings):
    localisation: MicroserviceSettings = MicroserviceSettings()
    esmfold: MicroserviceSettings = MicroserviceSettings()
    esmfold_light: MicroserviceSettings = MicroserviceSettings()
    rosettafold: MicroserviceSettings = MicroserviceSettings()
    gene_ontology: MicroserviceSettings = MicroserviceSettings()
    solubility: MicroserviceSettings = MicroserviceSettings()
    reinvent: MicroserviceSettings = MicroserviceSettings()
    protein_design: MicroserviceSettings = MicroserviceSettings()
    conformations: MicroserviceSettings = MicroserviceSettings()
    p2rank: MicroserviceSettings = MicroserviceSettings()
    msa_light: MsaLightMicroserviceSettings = MsaLightMicroserviceSettings()
    umol: MicroserviceSettings = MicroserviceSettings()
    diffdock: MicroserviceSettings = MicroserviceSettings()
    connection_string: str

    @classmethod
    def load(cls):
        return cls.parse_file(settings_path)

    def get_environment(self) -> Environment:
        key = 'NOLABS_ENVIRONMENT'
        if key not in os.environ:
            return Environment.LOCAL

        return Environment[key]
