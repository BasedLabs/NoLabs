import pandas as pd
from dockstream.core.result_parser import ResultParser

from dockstream.utils.enums.Gold_enums import GoldOutputEnum, GoldDockingConfigurationEnum


class GoldResultParser(ResultParser):
    """Class that loads, parses and analyzes the output of a "Gold" docking run, including poses and scores."""
    def __init__(self, ligands: list, fitness_function: str, response_value="fitness"):
        super().__init__(ligands=ligands)
        self._ROE = GoldOutputEnum()
        self._DE = GoldDockingConfigurationEnum()

        self._fitness_function = fitness_function
        self._response_value = response_value
        self._df_results = self._construct_dataframe()

    def _get_scoring_function_parameters(self):
        # get the appropriate name of the tag and whether mininmal or maximal values are best for
        # the specified scoring function
        if self._response_value == self._DE.GOLD_RESPONSE_VALUE_FITNESS:
            scoring_function_parameters = self._ROE.DICT_FITNESS[self._fitness_function]
        elif self._response_value == self._DE.GOLD_RESPONSE_VALUE_VALUE:
            scoring_function_parameters = self._ROE.DICT_VALUE[self._fitness_function]
        else:
            raise ValueError("Parameter response value must be either fitness or value.")
        self._logger.log(f"Set scoring_function_parameters to {scoring_function_parameters} for result parsing.",
                         self._LE.DEBUG)
        return scoring_function_parameters

    def _construct_dataframe(self) -> pd.DataFrame:
        scoring_function_parameters = self._get_scoring_function_parameters()

        def func_get_score(conformer):
            return float(conformer.GetProp(scoring_function_parameters[self._ROE.TAG]))

        return super()._construct_dataframe_with_funcobject(func_get_score)
