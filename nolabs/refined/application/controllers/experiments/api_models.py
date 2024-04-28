__all__ = [
    'ExperimentMetadataResponse',
    'ChangeExperimentNameRequest',
    'ExperimentMetadataResponse',
    'TimelineResponse'
]

import datetime
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class ExperimentMetadataResponse:
    id: UUID
    name: str
    date: datetime.datetime


@dataclass
class ChangeExperimentNameRequest:
    id: UUID
    name: str


@dataclass
class ExperimentMetadataResponse:
    id: UUID
    name: str
    date: datetime.datetime


@dataclass
class TimelineResponse:
    message: str
    error: str | None
    created_at: datetime.datetime
