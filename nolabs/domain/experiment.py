import dataclasses
import datetime
from typing import Dict

from pydantic.dataclasses import dataclass


@dataclasses.dataclass
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


@dataclasses.dataclass
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


@dataclasses.dataclass
@dataclass
class ExperimentMetadata:
    id: ExperimentId
    name: ExperimentName
    date: datetime.datetime
    properties: Dict
