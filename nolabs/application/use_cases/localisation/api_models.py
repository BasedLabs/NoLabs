from typing import List, Optional
from uuid import UUID

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
    protein_ids: List[UUID]
    result: List[JobResult]
    experiment_id: UUID


@dataclass
class SetupJobRequest:
    protein_ids: List[UUID]
    job_id: Optional[UUID] = None
    job_name: Optional[str] = None
    experiment_id: Optional[UUID] = None


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool
    result_valid: bool
