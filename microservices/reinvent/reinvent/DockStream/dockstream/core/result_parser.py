import abc
import pandas as pd
import warnings
from copy import deepcopy
from dockstream.core.ligand.ligand import Ligand

from dockstream.utils.dockstream_exceptions import ResultParsingFailed

from dockstream.loggers.docking_logger import DockingLogger
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.enums.docking_enum import ResultKeywordsEnum


class ResultParser(metaclass=abc.ABCMeta):
    """Base class implementing the interface result parsing classes."""

    def __init__(self, ligands):
        self._LE = LoggingConfigEnum()
        self._logger = DockingLogger()
        self._RK = ResultKeywordsEnum()

        self._ligands = ligands
        self._df_results = None

    def as_dataframe(self, aggregate=False):
        if aggregate:
            warnings.warn("For now, \"aggregate\" is not available.")
        if isinstance(aggregate, bool) and aggregate is False:
            return deepcopy(self._df_results)
        else:
            raise ResultParsingFailed("Parameter aggregate has an illegal value.")

    @staticmethod
    def _get_name(ligand: Ligand, conformer_index: int):
        """Function to make get either the name (for named molecules) or the identifier (plus the conformer) for the dataframe."""
        if ligand.get_name() is None:
            return ligand.get_identifier() + ':' + str(conformer_index)
        else:
            return ligand.get_name()

    def _construct_dataframe_with_funcobject(self, func_get_score) -> pd.DataFrame:
        data_buffer = []
        for ligand in self._ligands:
            best = True
            for conformer_index, conformer in enumerate(ligand.get_conformers()):
                name = self._get_name(ligand, conformer_index)
                row = [ligand.get_ligand_number(),
                       ligand.get_enumeration(),
                       conformer_index,
                       name,
                       func_get_score(conformer),
                       ligand.get_smile(),
                       best]
                best = False
                data_buffer.append(row)
        return pd.DataFrame(data_buffer, columns=[self._RK.DF_LIGAND_NUMBER,
                                                  self._RK.DF_LIGAND_ENUMERATION,
                                                  self._RK.DF_CONFORMER,
                                                  self._RK.DF_LIGAND_NAME,
                                                  self._RK.DF_SCORE,
                                                  self._RK.DF_SMILES,
                                                  self._RK.DF_LOWEST_CONFORMER])
