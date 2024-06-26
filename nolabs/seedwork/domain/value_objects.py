__all__ = [
    'ValueObject',
    'ValueObjectUUID',
    'ValueObjectString',
    'ValueObjectFloat',
    'ValueObjectBinary'
]

import uuid

from pydantic.dataclasses import dataclass


class ValueObject:
    """
    Base class for value objects
    """
    ...


@dataclass
class ValueObjectUUID(ValueObject):
    value: uuid.UUID

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, ValueObjectUUID):
            return False

        return self.value == other.value

    def __str__(self):
        return str(self.value)


@dataclass
class ValueObjectString(ValueObject):
    value: str

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, ValueObjectString):
            return False

        return self.value == other.value

    def __str__(self):
        return str(self.value)


@dataclass
class ValueObjectFloat(ValueObject):
    value: float

    def __float__(self):
        return self.value

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, ValueObjectFloat):
            return False

        return self.value == other.value

    def __str__(self):
        return str(self.value)


@dataclass
class ValueObjectBinary(ValueObject):
    value: bytes
