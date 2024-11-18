from typing import Literal

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    openai_api_key: str
    tavily_api_key: str


settings = Settings(_env_file=".env", _env_file_encoding="utf-8")
