import abc

from dockstream.utils.dockstream_exceptions import LigandPreparationFailed

from dockstream.loggers.ligand_preparation_logger import LigandPreparationLogger
from dockstream.loggers.blank_logger import BlankLogger
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.enums.transformations_enums import TransformationEnum


class Transformator(metaclass=abc.ABCMeta):
    """Base class implementing the interface for transformations applied to smiles before embedding."""

    def __init__(self, conf):
        self._LE = LoggingConfigEnum()
        self._TE = TransformationEnum()
        self._logger = LigandPreparationLogger()
        self._logger_blank = BlankLogger()
        self._conf = conf

        # extract type specification
        self._backend = conf[self._TE.TRANSFORMATION_BACKEND]
        if self._backend not in [self._TE.TRANSFORMATION_BACKEND_OPENEYE]:
            self._logger.log(f"Transformation backend {self._backend} is unknown.", self._LE.ERROR)
            raise LigandPreparationFailed(f"Transformation backend {self._backend} is unknown.")

        # extract backend specification
        self._type = conf[self._TE.TRANSFORMATION_TYPE]
        if self._type == self._TE.TRANSFORMATION_TYPE_SMIRKS:
            self._smirk = conf[self._TE.TRANSFORMATION_SMIRKS]
        else:
            self._logger.log(f"Transformation type {self._type} is unknown.", self._LE.ERROR)
            raise LigandPreparationFailed(f"Transformation type {self._type} is unknown.")

        # treat fail action specification
        self._fail_action = conf[self._TE.TRANSFORMATION_FAIL_ACTION]
        if self._fail_action not in [self._TE.TRANSFORMATION_FAIL_ACTION_KEEP,
                                     self._TE.TRANSFORMATION_FAIL_ACTION_DISCARD]:
            self._logger.log(f"Fail action {self._fail_action} is unknown.", self._LE.ERROR)
            raise LigandPreparationFailed(f"Fail action {self._fail_action} is unknown.")

    def transform(self, ligands: list) -> list:
        raise NotImplementedError
