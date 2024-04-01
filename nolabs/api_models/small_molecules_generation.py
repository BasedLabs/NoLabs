import datetime

from pydantic import dataclasses

@dataclasses.dataclass
class JobResponse:
    id: str
    name: str
    created_at: datetime.datetime
    running: bool
    learning_completed: bool


@dataclasses.dataclass
class LogsResponse:
    output: str
    docking_output: str
    errors: str


@dataclasses.dataclass
class ParamsResponse:
    pdb_file: str
    pdb_file_name: str
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    batch_size: float
    minscore: float
    epochs: float


@dataclasses.dataclass
class SmilesResponse:
    smiles: str
    drug_likeness: float
    score: float
    created_at: datetime.datetime