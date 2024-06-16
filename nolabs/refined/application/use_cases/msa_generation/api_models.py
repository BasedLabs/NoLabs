__all__ = [
    'JobResponse',
    'SetupJobRequest',
    'RunJobRequest',
    'GetJobStatusResponse'
]

from typing import Optional
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    experiment_id: UUID
    protein_id: UUID
    result: str | None


@dataclass
class SetupJobRequest:
    experiment_id: UUID
    protein_id: UUID

    job_id: Optional[UUID] = None
    job_name: Optional[str] = None


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool


