__all__ = [
    'JobResult',
    'JobResponse',
    'RunJobRequest',
    'GetJobStatusResponse'
]


from enum import Enum
from typing import List, Any, Optional
from uuid import UUID

from pydantic import BaseModel, model_validator
from pydantic.dataclasses import dataclass


@dataclass
class JobResult:
    ligand_id: UUID
    sdf_content: str
    minimized_affinity: Optional[float]
    scored_affinity: Optional[float]
    confidence: Optional[float]


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    samples_per_complex: int
    protein_id: UUID
    ligand_ids: List[UUID]
    result: List[JobResult]


@dataclass
class SetupJobRequest:
    job_id: Optional[UUID]
    job_name: Optional[str]
    experiment_id: UUID
    protein_id: UUID
    ligand_ids: List[UUID]
    samples_per_complex: Optional[int] = 40


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool
