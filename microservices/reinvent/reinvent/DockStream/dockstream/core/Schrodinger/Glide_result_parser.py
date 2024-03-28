import pandas as pd
from dockstream.core.result_parser import ResultParser

from dockstream.utils.enums.Schrodinger_enums import SchrodingerOutputEnum


class GlideResultParser(ResultParser):
    """Class that loads, parses and analyzes the output of a "Glide" docking run, including poses and scores."""
    def __init__(self, ligands: list):
        super().__init__(ligands=ligands)

        self._df_results = self._construct_dataframe()

    def _construct_dataframe(self) -> pd.DataFrame:
        def func_get_score(conformer):
            _ROE = SchrodingerOutputEnum()
            return float(conformer.GetProp(_ROE.GLIDE_DOCKING_SCORE))

        return super()._construct_dataframe_with_funcobject(func_get_score)
