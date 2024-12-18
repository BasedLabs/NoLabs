from enum import Enum
from pathlib import Path
from typing import Literal, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

directory = Path(__file__).resolve().parent


class Environment(str, Enum):
    local = 'local'
    test = 'test'
    production = 'production'


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=[directory / ".env"],
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
    redis_url: str
    workflow_queue: str = "workflow"
    celery_worker_pool: str = "prefork"
    celery_worker_state_db: str = "/tmp/celery-state.db"
    connection_string: str
    enable_biobuddy: bool
    reinvent_directory: Path
    uvicorn_host: str = "0.0.0.0"
    adaptyv_bio_api_token: Optional[str] = None
    adaptyv_bio_api_base: Optional[str] = None
    uvicorn_port: int = 8000
    celery_worker_concurrency: int = 20
    environment: Literal["local", "test", "production"] = "local"
    logging_level: Literal["INFO", "WARNING", "ERROR"] = "INFO"
    structured_logging: bool = False


_settings: Optional[Settings] = None


class SettingsProxy:
    def __getattr__(self, name):
        global _settings
        return getattr(_settings, name)

    def __call__(self):
        global _settings
        return self._settings


settings: Settings = SettingsProxy()  # type: ignore


def initialize_settings():
    global _settings

    _settings = Settings()
