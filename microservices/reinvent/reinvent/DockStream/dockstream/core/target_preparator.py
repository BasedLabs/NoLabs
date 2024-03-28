import abc

from dockstream.containers.target_preparation_container import TargetPreparationContainer
from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum

from dockstream.loggers.target_preparation_logger import TargetPreparationLogger
from dockstream.loggers.blank_logger import BlankLogger

from dockstream.utils.enums.logging_enums import LoggingConfigEnum


class TargetPreparator(metaclass=abc.ABCMeta):
    """Virtual base class implementing the interface for all specific target preparators and the general preparation
       of the docking target."""

    def __init__(self, conf: TargetPreparationContainer, run_number=0):
        self._TE = TargetPreparationEnum()
        self._TL = LoggingConfigEnum()
        self._logger = TargetPreparationLogger()
        self._logger_blank = BlankLogger()
        self._config = conf
        self._target = None

        # store the specific parameters for this very run for easy access later; the others are ignored
        self._run_parameters = self._config[self._TE.TARGETPREP][self._TE.RUNS][run_number]

    def get_target(self):
        return self._target

    def specify_cavity(self):
        raise NotImplementedError("This method needs to be overwritten by child classes.")

    def write_target(self, path):
        raise NotImplementedError("This method needs to be overwritten by child classes.")
