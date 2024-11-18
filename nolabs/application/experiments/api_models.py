__all__ = [
    "ExperimentMetadataResponse",
    "UpdateExperimentRequest",
    "TimelineResponse",
]

import datetime
from typing import Optional
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
class TimelineResponse:
    created_at: datetime.datetime
    message: Optional[str] = None
    error: Optional[str] = None
