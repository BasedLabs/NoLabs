from typing import Literal

from pydantic import (
    computed_field,
)
from pydantic_settings import BaseSettings, SettingsConfigDict

from infrastructure.environment import Environment


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=["infrastructure/.env"],
        case_sensitive=False,
        env_ignore_empty=True,
        extra="ignore",
        env_prefix='NOLABS_'
    )
    localisation_host: str
    biobuddy_host: str
    external_query_host: str
    esmfold_host: str
    esmfold_light_host: str
    rosettafold_host: str
    gene_ontology_host: str
    solubility_host: str
    reinvent_host: str
    protein_design_host: str
    conformations_host: str
    p2rank_host: str
    msa_light_host: str
    msa_light_server_url: str
    umol_host: str
    diffdock_host: str
    redis_host: str
    redis_port: str
    celery_broker: str
    celery_backend: str
    domain: str = "localhost"
    connection_string: str
    socketio_broker: str
    environment: Literal["local", "test", "production"] = "local"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def server_host(self) -> str:
        if self.environment == "local":
            return f"http://{self.domain}"
        return f"https://{self.domain}"

    def get_environment(self) -> Environment:
        return Environment[self.environment.upper()]


settings = Settings()  # type: ignore
