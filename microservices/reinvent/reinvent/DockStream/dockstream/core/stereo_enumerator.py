from abc import ABC

from pydantic import BaseModel, PrivateAttr

from dockstream.loggers.ligand_preparation_logger import LigandPreparationLogger
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.enums.stereo_enumeration_enums import StereoEnumerationEnum

_LE = LoggingConfigEnum()
_SE = StereoEnumerationEnum()


class StereoEnumerator(ABC, BaseModel):
    _logger = PrivateAttr()

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)
        self._logger = LigandPreparationLogger()

    def enumerate(self, ligands: list) -> list:
        raise NotImplementedError
