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
