import configparser
import os


class Settings:
    def __init__(self):
        self._config = configparser.ConfigParser()
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self._config.read(os.path.join(dir_path, 'settings.ini'))

    @property
    def model_name(self) -> str:
        return self._config.get('model', 'model_name')
