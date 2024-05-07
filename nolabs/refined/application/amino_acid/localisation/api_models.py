from typing import List, Optional, Any
from uuid import UUID

from pydantic import model_validator, BaseModel
from pydantic.dataclasses import dataclass


@dataclass
class JobResult:
    protein_id: UUID
    cytosolic: float
    mitochondrial: float
    nuclear: float
    other: float
    extracellular: float


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    proteins: List[UUID]
    result: List[JobResult]


@dataclass
class SetupJobRequest:
    job_id: Optional[UUID]
    job_name: Optional[str]
    experiment_id: UUID
    proteins: List[UUID]


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool
