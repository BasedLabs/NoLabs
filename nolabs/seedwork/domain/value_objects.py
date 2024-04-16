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

    @classmethod
    def validate(cls, value):
        try:
            uuid.UUID(str(value))
            return True
        except ValueError:
            return False


@dataclass
class ValueObject:
    """
    Base class for value objects
    """
    ...


TReadableContent = TypeVar('TReadableContent')


@dataclass
class ReadableValueObject(ValueObject, Generic[TReadableContent], ABC):
    @abstractmethod
    def read(self) -> TReadableContent:
        ...

    @abstractmethod
    def write(self, content: TReadableContent):
        ...
