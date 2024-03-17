from __future__ import annotations

import datetime
from typing import List

import pydantic


@pydantic.dataclasses.dataclass
class RunReinventResponse:
    pdb_content: str | None
    errors: List[str]


@pydantic.dataclasses.dataclass
class RunFineTuningJobRequest:
    name: str
    epochs: int = 50


@pydantic.dataclasses.dataclass
class FineTuningJobResponse:
    id: str
    name: str
    running: bool
    started_at: datetime.datetime
    finished_at: datetime.datetime
    progress: float
    pdb_content: str
    pdb_filename: str
    errors: str
    epochs: int
