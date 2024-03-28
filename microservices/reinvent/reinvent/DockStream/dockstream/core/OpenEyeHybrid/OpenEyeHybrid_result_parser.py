import pandas as pd
from dockstream.core.result_parser import ResultParser

from dockstream.utils.enums.OE_Hybrid_enums import OpenEyeHybridOutputKeywordsEnum


class OpenEyeHybridResultParser(ResultParser):
    """Loads, parses and analyzes the output of an "OpenEye Hybrid" docking run, including poses and score."""
    def __init__(self, ligands: list):
        super().__init__(ligands=ligands)
        self._OE = OpenEyeHybridOutputKeywordsEnum()
        self._df_results = self._construct_dataframe()

    def _construct_dataframe(self) -> pd.DataFrame:
        def func_get_score(conformer):
            return float(conformer.GetProp(self._OE.SCORE))

        return super()._construct_dataframe_with_funcobject(func_get_score)
