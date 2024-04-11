from enum import IntEnum


class ErrorCode(IntEnum):
    PDBQT_NOT_PROVIDED = 1
    PDB_NOT_PROVIDED = 2
    LEARNING_IS_NOT_COMPLETED = 3
    NUMBER_OF_MOLECULES_REQUIRED = 4


class ReinventException(Exception):
    def __init__(self, code: ErrorCode):
        self.code: ErrorCode = code
