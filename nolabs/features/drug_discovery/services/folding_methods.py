from enum import Enum

from nolabs.exceptions import ErrorCodes, NoLabsException


class FoldingMethods(str, Enum):
    esmfold = 'esmfold'
    esmfold_light = 'esmfold_light'
    rosettafold = 'rosettafold'

    @staticmethod
    def ensure_contains(s: str):
        """
        :raises NoLabsException: if parameter is not a valid FoldingMethods enum
        """
        if s not in [e for e in FoldingMethods]:
            raise NoLabsException([s], ErrorCodes.folding_method_unknown)
