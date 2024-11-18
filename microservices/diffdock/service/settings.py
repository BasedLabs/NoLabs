from pathlib import Path
from typing import Literal

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    redis_url: str
    celery_worker_queue: str
    celery_worker_concurrency: int
    fastapi_host: str
    fastapi_port: int
    weights_path: Path
    home: Path
    use_max_power: bool
    omp_num_threads: int
    mkl_num_threads: int
    logging_level: Literal["DEBUG", "WARNING", "ERROR", "INFO"] = "INFO"
    environment: Literal["local", "test", "production"] = "local"


settings = Settings(_env_file=".env", _env_file_encoding="utf-8")
