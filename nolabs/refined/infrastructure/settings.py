__all__ = [
    'Settings'
]

import os
from enum import Enum

from pydantic import Field, BaseModel
from pydantic_settings import BaseSettings, SettingsConfigDict


class Environment(str, Enum):
    LOCAL = 'local'
    TEST = 'test'


class MicroserviceSettings(BaseModel):
    microservice: str = ''


settings_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                             f'appsettings.{os.environ["NOLABS_ENVIRONMENT"]}.json')


class Settings(BaseSettings):
    localisation: MicroserviceSettings = MicroserviceSettings()
    connection_string: str

    @classmethod
    def load(cls):
        return cls.parse_file(settings_path)

    def get_environment(self) -> Environment:
        key = 'NOLABS_ENVIRONMENT'
        if key not in os.environ:
            return Environment.LOCAL

        return Environment[key]
