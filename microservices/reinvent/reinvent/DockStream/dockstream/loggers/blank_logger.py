import logging

from dockstream.loggers.base_logger import BaseLogger


class BlankLogger(BaseLogger):
    """This logger serves as a "verbatim" interface."""
    def __init__(self):
        super().__init__()

    def _initialize_logger(self):
        logger = logging.getLogger(self._LE.LOGGER_BLANK)
        return logger
