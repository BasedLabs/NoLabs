from __future__ import annotations

from enum import IntEnum


class ErrorCodes(IntEnum):
    MUST_BE_ABS_PATH = 1
    ALREADY_EXISTS = 2
    OBJECT_DOES_NOT_EXIST = 3
    PDB_FILE_DOES_NOT_EXIST = 4


class LeafException(BaseException):
    def __init__(self, code: ErrorCodes):
        self.message = ''
        self.code = code

    def with_additional_message(self, s: str) -> LeafException:
        self.message = s
        return self
