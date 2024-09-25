__all__ = ["JobResult", "JobResponse", "RunJobRequest"]


from typing import List, Optional
from uuid import UUID

from pydantic import BaseModel, model_validator
from pydantic.dataclasses import dataclass


@dataclass
class JobResult:
    complex_id: UUID
    sdf_content: str
    minimized_affinity: Optional[float] = None
    scored_affinity: Optional[float] = None
    confidence: Optional[float] = None


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    samples_per_complex: int
    protein_id: UUID
    ligand_id: UUID
    result: List[JobResult]


@dataclass
class SetupJobRequest:
    experiment_id: UUID
    protein_id: UUID
    ligand_id: UUID
    samples_per_complex: Optional[int] = 2
    job_id: Optional[UUID] = None
    job_name: Optional[str] = None


@dataclass
class RunJobRequest:
    job_id: UUID
