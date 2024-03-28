import pandas as pd
from dockstream.core.result_parser import ResultParser

from dockstream.utils.enums.OpenEye_enums import OpenEyeResultKeywordsEnum


class OpenEyeResultParser(ResultParser):
    """Class that loads, parses and analyzes the output of an "OpenEye" docking run, including poses and scores."""
    def __init__(self, ligands: list):
        super().__init__(ligands=ligands)
        self._RK = OpenEyeResultKeywordsEnum()

        self._df_results = self._construct_dataframe()

    def _construct_dataframe(self) -> pd.DataFrame:
        def func_get_score(conformer):
            return float(conformer.GetEnergy())

        return super()._construct_dataframe_with_funcobject(func_get_score)
