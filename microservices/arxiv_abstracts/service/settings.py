from typing import Literal, Optional

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    redis_url: str
    celery_enable_utc: bool
    celery_worker_concurrency: int
    celery_worker_queue: str
    fastapi_host: str
    fastapi_port: int
    chroma_db_path: str
    logging_level: Literal["DEBUG", "WARNING", "ERROR", "INFO"] = "INFO"
    environment: Literal["local", "test", "production"] = "local"
    openai_api_key: Optional[str] = None


settings = Settings(_env_file=".env", _env_file_encoding="utf-8")
