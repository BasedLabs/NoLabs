from enum import IntEnum


class ErrorCode(IntEnum):
    PDBQT_NOT_PROVIDED = 1


class ReinventException(Exception):
    def __init__(self, code: ErrorCode):
        self.code: ErrorCode = code
