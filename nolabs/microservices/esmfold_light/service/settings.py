from typing import Literal

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    celery_broker_url: str
    celery_backend_url: str
    celery_enable_utc: bool
    fastapi_host: str
    fastapi_port: int
    logging_level: Literal['debug', 'info'] = 'info'
    mode: Literal['fastapi', 'celery'] = 'celery'


settings = Settings(_env_file='.env', _env_file_encoding='utf-8')
