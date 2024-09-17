from pathlib import Path
from typing import Literal

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    celery_broker_url: str
    celery_backend_url: str
    celery_enable_utc: bool
    celery_worker_concurrency: int
    fastapi_host: str
    fastapi_port: int
    configs_directory: Path
    logging_level: Literal["DEBUG", "WARNING", "ERROR", "INFO"] = "INFO"
    mode: Literal["fastapi", "celery"] = "celery"
    environment: Literal["local", "test", "production"] = "local"


settings = Settings(_env_file=".env", _env_file_encoding="utf-8")
