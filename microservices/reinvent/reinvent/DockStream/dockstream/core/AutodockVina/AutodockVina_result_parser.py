import pandas as pd
from dockstream.core.result_parser import ResultParser

from dockstream.utils.enums.AutodockVina_enums import AutodockResultKeywordsEnum


class AutodockResultParser(ResultParser):
    """Class that loads, parses and analyzes the output of an "AutoDock Vina" docking run, including poses and scores."""
    def __init__(self, ligands):
        super().__init__(ligands=ligands)
        self._RK = AutodockResultKeywordsEnum()

        self._df_results = self._construct_dataframe()

    def _construct_dataframe(self) -> pd.DataFrame:
        def func_get_score(conformer):
            return float(conformer.GetProp(self._RK.SDF_TAG_SCORE))

        return super()._construct_dataframe_with_funcobject(func_get_score)
