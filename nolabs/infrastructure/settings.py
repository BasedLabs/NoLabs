from pathlib import Path
from typing import Literal, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

directory = Path(__file__).resolve().parent


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
    celery_broker: str
    celery_backend: str
    workflow_queue: str = "workflow"
    celery_worker_pool: str = "prefork"
    celery_worker_state_db: str = "/tmp/celery-state.db"
    connection_string: str
    socketio_broker: str
    enable_structured_logging: bool
    reinvent_directory: Path
    blast_email: str
    uvicorn_host: str = "0.0.0.0"
    uvicorn_port: int = 8000
    celery_worker_concurrency: int = 20
    orphaned_tasks_check_interval: int = 40
    environment: Literal["local", "test", "production"] = "local"
    logging_level: Literal["INFO", "WARNING", "ERROR"] = "INFO"


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
