__all__ = [
    'ExperimentMetadataResponse',
    'UpdateExperimentRequest',
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
class UpdateExperimentRequest:
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
