from enum import IntEnum


class ErrorCode(IntEnum):
    PDB_NOT_PROVIDED = 1


class ReinventException(Exception, IntEnum):
    def __init__(self, code: ErrorCode):
        self.code: ErrorCode = code
