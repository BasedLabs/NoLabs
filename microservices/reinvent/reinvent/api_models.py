from __future__ import annotations

import datetime
from typing import List

import pydantic


@pydantic.dataclasses.dataclass
class ParamsRequest:
    config_id: str

    # size of the search box. Recommended not more than 30 armstrongs
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    batch_size: int
    minscore: float
    name: str
    epochs: int = 50



@pydantic.dataclasses.dataclass
class ConfigurationResponse:
    id: str
    name: str
    created_at: datetime.datetime
    running: bool
    sampling_allowed: bool


@pydantic.dataclasses.dataclass
class LogsResponse:
    output: str
    docking_output: str
    errors: str


@pydantic.dataclasses.dataclass
class ParamsResponse:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    batch_size: float
    minscore: float
    epochs: float


@pydantic.dataclasses.dataclass
class Smiles:
    smiles: str
    drugLikeness: float
    score: float
    stage: str


@pydantic.dataclasses.dataclass
class SmilesResponse:
    smiles: List[Smiles]


@pydantic.dataclasses.dataclass
class SamplingSizeRequest:
    number_of_molecules_to_design: int
