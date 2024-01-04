from enum import Enum
from typing import List


class ExceptionCodes(Enum):
    fix_pdb_error = 1


class NoLabsException(Exception):
    def __init__(self, message: List[str] | str, error_code: ExceptionCodes):
        self.message = '|'.join(message) if isinstance(message, str) else message
        self.error_code = error_code
