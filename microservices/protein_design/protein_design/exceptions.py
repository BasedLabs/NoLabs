from enum import Enum
from typing import List


class ExceptionCodes(Enum):
    rfdiffusion_folder_does_not_exists = 1


class ProteinDesignException(Exception):
    def __init__(self, message: List[str] | str, error_code: ExceptionCodes):
        self.message = '|'.join(message) if isinstance(message, str) else message
        self.error_code = error_code