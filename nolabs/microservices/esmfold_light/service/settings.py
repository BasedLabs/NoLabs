from typing import Literal

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    celery_broker_url: str
    celery_backend_url: str
    celery_enable_utc: bool
    fastapi_host: str
    fastapi_port: int
    application_loglevel: Literal['debug', 'info'] = 'info'
    mode: Literal['fastapi', 'celery'] = 'celery'


settings = Settings(_env_file='.env.development', _env_file_encoding='utf-8')
