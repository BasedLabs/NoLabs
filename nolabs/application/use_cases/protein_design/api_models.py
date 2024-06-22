from __future__ import annotations
from typing import List, Optional
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    experiment_id: UUID
    protein_id: UUID

    binder_ids: List[UUID]

    contig: Optional[str] = None
    number_of_designs: Optional[int] = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None


@dataclass
class SetupJobRequest:
    experiment_id: UUID
    protein_id: UUID

    contig: str
    number_of_designs: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None
    job_id: Optional[UUID] = None
    job_name: Optional[str] = None


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool
    result_valid: bool
