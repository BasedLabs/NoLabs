import datetime
from typing import Optional
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class SmilesResponse:
    smiles: str
    drug_likeness: float
    score: float
    created_at: datetime.datetime
    stage: str


@dataclass
class LogsResponse:
    output: str
    docking_output: str
    errors: str


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    experiment_id: UUID
    protein_id: UUID

    # Learning
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    batch_size: int
    minscore: float
    epochs: int

    # Sampling
    sampling_size: int


@dataclass
class StartSamplingRequest:
    sampling_size: int = 5


@dataclass
class SetupJobRequest:
    protein_id: UUID

    # Learning
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    batch_size: int = 50
    minscore: float = 0.4
    epochs: int = 128

    job_id: Optional[UUID] = None
    job_name: Optional[str] = None

    # Sampling
    sampling_size: int = 5

    experiment_id: Optional[UUID] = None


@dataclass
class GetJobStatusResponse:
    running: bool
    sampling_allowed: bool
    result_valid: bool
