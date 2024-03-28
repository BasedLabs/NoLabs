from dockstream.loggers.ligand_preparation_logger import LigandPreparationLogger
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.core.OpenEye.OpenEye_transformator import OpenEyeTransformator
from dockstream.utils.enums.transformations_enums import TransformationEnum


class TransformatorFactory:
    """Returns a list of transformators."""

    def __init__(self, conf):
        self._TE = TransformationEnum()
        self._LE = LoggingConfigEnum()
        self._logger = LigandPreparationLogger()
        self._conf = conf

    def get_transformators(self) -> list:
        transformators = []
        for curTransConf in self._conf[self._TE.TRANSFORMATIONS]:
            if curTransConf[self._TE.TRANSFORMATION_BACKEND] == self._TE.TRANSFORMATION_BACKEND_OPENEYE:
                transformators.append(OpenEyeTransformator(curTransConf))
            else:
                self._logger.log(f"", self._LE.DEBUG)

        return transformators
