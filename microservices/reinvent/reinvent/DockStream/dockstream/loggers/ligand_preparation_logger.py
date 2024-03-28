import logging

from dockstream.loggers.base_logger import BaseLogger


class LigandPreparationLogger(BaseLogger):
    def __init__(self):
        super().__init__()

    def _initialize_logger(self):
        logger = logging.getLogger(self._LE.LOGGER_LIGAND_PREPARATION)
        return logger
