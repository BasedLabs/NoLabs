import datetime
from typing import Optional

from fastapi import UploadFile
from pydantic import dataclasses


@dataclasses.dataclass
class LogsResponse:
    output: str
    docking_output: str
    errors: str


@dataclasses.dataclass
class ExperimentPropertiesResponse:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    batch_size: float
    minscore: float
    epochs: float
    pdb_file: Optional[str]
    pdb_file_name: Optional[str]


@dataclasses.dataclass
class GetExperimentStatusResponse:
    running: bool
    sampling_allowed: bool


@dataclasses.dataclass
class GetExperimentResponse:
    experiment_id: str
    experiment_name: str
    created_at: datetime.datetime
    status: GetExperimentStatusResponse
    properties: ExperimentPropertiesResponse


@dataclasses.dataclass
class ExperimentPropertiesRequest:
    pdb_file: UploadFile

    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    batch_size: int
    minscore: float
    epochs: int = 50


@dataclasses.dataclass
class SmilesResponse:
    smiles: str
    drug_likeness: float
    score: float
    created_at: datetime.datetime
    stage: str

@dataclasses.dataclass
class SamplingSizeRequest:
    number: int