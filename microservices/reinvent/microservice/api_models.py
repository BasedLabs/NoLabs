from __future__ import annotations

import datetime

import pydantic


@pydantic.dataclasses.dataclass
class RunFineTuningJobRequest:
    # size of the search box. Recommended not more than 30 armstrongs
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    epochs: int = 50


@pydantic.dataclasses.dataclass
class FineTuningJobResponse:
    id: str
    running: bool
    started_at: datetime.datetime
    progress: float
    pdbqt_content: str
    pdbqt_filename: str
    errors: str
    epochs: int


@pydantic.dataclasses.dataclass
class RunInferenceRequest:
    job_id: str


@pydantic.dataclasses.dataclass
class RunInferenceResponse:
    id: str
    job_name: str
    smiles: str | None
    errors: str | None
