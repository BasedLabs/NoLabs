from typing import Literal

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    celery_broker_url: str
    celery_backend_url: str
    celery_enable_utc: bool
    celery_worker_concurrency: int
    celery_worker_queue: str
    fastapi_host: str
    fastapi_port: int
    rfdiffusion_path: str
    logging_level: Literal["DEBUG", "WARNING", "ERROR", "INFO"] = "INFO"
    environment: Literal["local", "test", "production"] = "local"


settings = Settings(_env_file=".env", _env_file_encoding="utf-8")
