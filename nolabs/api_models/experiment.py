import datetime

import pydantic


@pydantic.dataclasses.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str


@pydantic.dataclasses.dataclass
class ExperimentMetadataResponse:
    id: str
    name: str
    date: datetime.datetime


@pydantic.dataclasses.dataclass
class TimelineResponse:
    message: str
    error: str | None
    created_at: datetime.datetime
