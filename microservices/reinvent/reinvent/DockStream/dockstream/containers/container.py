import abc
import json
import os


class ConfigurationContainer(object, metaclass=abc.ABCMeta):
    """Class that takes a lot of arguments in the form of a JSON configuration (as dictionary, string or file path)
       and, optionally, performs a JSON Schema validation."""

    @abc.abstractmethod
    def __init__(self, conf):
        # get instance of configuration enum and load configuration
        # parameter "config" can be a string, a path or a dictionary (as long as it holds valid JSON input)
        if isinstance(conf, str):
            if os.path.isfile(conf):
                with open(conf) as file:
                    conf = file.read().replace("\r", "").replace("\n", "")
            conf = json.loads(conf)
        self._conf = conf

    def get_as_dict(self):
        return self._conf

    def get(self, key, default=None):
        return self._conf.get(key, default)

    def __getitem__(self, item):
        return self.get_as_dict()[item]

    def get_as_string(self):
        return json.dumps(self._conf)

    def validate(self):
        raise NotImplementedError("This functions needs to be implemented by child classes.")
