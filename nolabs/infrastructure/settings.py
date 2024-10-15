from pathlib import Path
from typing import Literal

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=[".env"],
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
    celery_worker_pool: str = 'prefork'
    connection_string: str
    socketio_broker: str
    enable_structured_logging: bool
    reinvent_directory: Path
    blast_email: str
    workflow_version: int
    celery_worker_concurrency: int = 1
    mode: Literal["united", "fastapi", "workflow"] = "fastapi"
    environment: Literal["local", "test", "production"] = "local"
    logging_level: Literal["INFO", "WARNING", "ERROR"] = "INFO"


settings = Settings()
