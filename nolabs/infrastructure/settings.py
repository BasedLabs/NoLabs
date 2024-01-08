import configparser
import os.path
from dataclasses import dataclass

__all__ = ['Settings', ]

from nolabs.infrastructure.environment import is_dev

import configparser
import os


class Settings:
    def __init__(self):
        self._config = configparser.ConfigParser()
        dir_path = os.path.dirname(os.path.realpath(__file__))

        settings_path = 'settings-dev.ini' if is_dev() else 'settings.ini'

        self._config.read(os.path.join(dir_path, settings_path))

    @property
    def conformations_simulations_file_name(self) -> str:
        return self._config.get('conformations', 'simulations_file_name')

    @property
    def conformations_experiments_folder(self) -> str:
        import nolabs
        exp_folder = self._config.get('conformations', 'conformations_experiments_folder')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def conformations_metadata_file_name(self) -> str:
        return self._config.get('conformations', 'conformations_metadata_file_name')
