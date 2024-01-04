import configparser
import os.path
from dataclasses import dataclass

__all__ = ['get_config', ]

from nolabs.infrastructure.environment import is_dev


@dataclass
class MicroservicesConfiguration:
    conformations_api_url: str

@dataclass
class Configuration:
    microservices_configuration: MicroservicesConfiguration


def get_config() -> Configuration:
    settings_file_name = ''

    if is_dev():
        settings_file_name = 'settings-dev.ini'

    assert settings_file_name

    configuration_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), settings_file_name)
    parser = configparser.ConfigParser()
    parser.read(configuration_path)

    microservices_configuration = MicroservicesConfiguration(parser['microservice'].get('conformations_api_url'))

    return Configuration(microservices_configuration)
