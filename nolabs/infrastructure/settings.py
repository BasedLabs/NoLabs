import configparser
import os.path
from dataclasses import dataclass

__all__ = ['get_config', ]

from nolabs.infrastructure.environment import is_dev


import configparser
import os


class Settings:
    def __init__(self):
        self._config = configparser.ConfigParser()
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self._config.read(os.path.join(dir_path, 'settings.ini'))

    @property
    def model_weights_name(self) -> str:
        return self._config.get('model', 'model_weights_name')
