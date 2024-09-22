from pathlib import Path
from typing import Literal

from pydantic import computed_field
from pydantic_settings import BaseSettings, SettingsConfigDict

from nolabs.infrastructure.environment import Environment


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=["infrastructure/.env"],
        case_sensitive=False,
        env_ignore_empty=True,
        extra="ignore",
        env_prefix="NOLABS_",
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
    connection_string: str
    socketio_broker: str
    enable_structured_logging: bool
    home: Path
    reinvent_directory: Path
    blast_email: str
    environment: Literal["local", "test", "production"] = "local"
    logging_level: Literal["INFO", "WARNING", "ERROR"] = "INFO"

    def get_environment(self) -> Environment:
        return Environment[self.environment.upper()]


settings = Settings()  # type: ignore

if not Path(settings.home).exists():
    raise RuntimeError(f"{settings.home} does not exist. Create it.")
