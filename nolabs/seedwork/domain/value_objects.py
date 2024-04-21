import uuid
from abc import abstractmethod, ABC
from typing import TypeVar, Generic

from pydantic import BaseModel, ConfigDict
from pydantic.dataclasses import dataclass


class GenericUUID(uuid.UUID):
    @classmethod
    def next_id(cls):
        return cls(int=uuid.uuid4().int)

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    def validate(cls, value):
        try:
            uuid.UUID(str(value))
            return True
        except ValueError:
            return False


class ValueObject:
    """
    Base class for value objects
    """
    ...