__all__ = [
    'Settings'
]

import os

from pydantic import Field, BaseModel
from pydantic_settings import BaseSettings, SettingsConfigDict


class MicroserviceSettings(BaseModel):
    microservice: str = Field()


class Settings(BaseSettings):
    model_config = SettingsConfigDict(json_file=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                             f'appsettings.{os.environ["NOLABS_ENVIRONMENT"]}.json'))

    localisation: MicroserviceSettings = Field()
