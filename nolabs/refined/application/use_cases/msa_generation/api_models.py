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
    msa: str | None


@dataclass
class SetupJobRequest:
    job_id: Optional[UUID]
    job_name: Optional[str]
    experiment_id: UUID
    protein_id: UUID


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool


