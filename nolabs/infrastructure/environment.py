from enum import Enum


class Environment(str, Enum):
    LOCAL = 'local'
    TEST = 'test'
    UNITTESTS = 'unit'

