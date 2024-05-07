from __future__ import annotations
from typing import List, Optional
from uuid import UUID

from fastapi import UploadFile
from pydantic.dataclasses import dataclass


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    experiment_id: UUID
    protein_id: UUID

    binder_ids: List[UUID]

    number_of_designs: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None


@dataclass
class SetupJobRequest:
    job_id: Optional[UUID]
    job_name: Optional[str]
    experiment_id: UUID
    protein_id: UUID

    binder_ids: List[UUID]

    number_of_designs: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool