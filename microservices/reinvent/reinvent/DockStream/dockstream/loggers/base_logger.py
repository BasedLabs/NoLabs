from abc import ABC, abstractmethod

from dockstream.utils.enums.logging_enums import LoggingConfigEnum


class BaseLogger(ABC):
    def __init__(self):
        self._LE = LoggingConfigEnum()
        self._logger = self._initialize_logger()

    def log(self, message: str, level: str):
        if level == self._LE.DEBUG:
            self._logger.debug(message)
        elif level == self._LE.INFO:
            self._logger.info(message)
        elif level == self._LE.WARNING:
            self._logger.warning(message)
        elif level == self._LE.ERROR:
            self._logger.error(message)
        elif level == self._LE.EXCEPTION:
            self._logger.exception(message)
        else:
            raise ValueError("Logger level not supported.")

    @abstractmethod
    def _initialize_logger(self):
        raise NotImplementedError("Overwrite this method in child classes.")
