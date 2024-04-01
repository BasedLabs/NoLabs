import datetime
from typing import TypeVar

from pydantic.dataclasses import dataclass

ExperimentPropertiesT = TypeVar('ExperimentPropertiesT')

@dataclass
class ExperimentId:
    value: str

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise ValueError()

        return self.value == other.value

    def __str__(self):
        return self.value

    def __hash__(self):
        return self.value.__hash__()


@dataclass
class ExperimentName:
    value: str

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise ValueError()

        return self.value == other.value

    def __str__(self):
        return self.value

    def __hash__(self):
        return self.value.__hash__()


@dataclass
class ExperimentMetadata:
    id: ExperimentId
    name: ExperimentName
    date: datetime.datetime
