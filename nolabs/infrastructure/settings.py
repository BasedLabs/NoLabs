from enum import Enum
from pathlib import Path
from typing import Literal

from pydantic_settings import BaseSettings, SettingsConfigDict

from nolabs.infrastructure.environment import Environment


class RunModeEnum(str, Enum):
    united = "united"


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=["infrastructure/.env"],
        case_sensitive=False,
        env_ignore_empty=True,
        extra="ignore",
        env_prefix="NOLABS_",
    )
    biobuddy_host: str
    external_query_host: str
    esmfold_light_host: str
    reinvent_host: str
    diffdock_host: str
    celery_broker: str
    celery_backend: str
    workflow_queue: str = 'workflow'
    connection_string: str
    socketio_broker: str
    enable_structured_logging: bool
    reinvent_directory: Path
    blast_email: str
    workflow_version: int
    mode: Literal["united", "fastapi", "workflow"] = "fastapi"
    environment: Literal["local", "test", "production"] = "local"
    logging_level: Literal["INFO", "WARNING", "ERROR"] = "INFO"

    def get_environment(self) -> Environment:
        return Environment[self.environment.upper()]


settings = Settings()  # type: ignore
