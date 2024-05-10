import datetime

from pydantic.dataclasses import dataclass


@dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str

@dataclass
class ExperimentMetadataRequest:
    id: str

@dataclass
class ExperimentMetadataResponse:
    id: str
    name: str
    date: datetime.datetime


@dataclass
class TimelineResponse:
    message: str
    error: str | None
    created_at: datetime.datetime
