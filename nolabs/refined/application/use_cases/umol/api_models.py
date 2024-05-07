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
    predicted_pdb: bytes
    predicted_sdf: bytes
    plddt_array: List[int]


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    pocket_ids: List[int]
    protein_id: UUID
    ligand_id: UUID
    result: JobResult


@dataclass
class SetupJobRequest:
    job_id: Optional[UUID]
    job_name: Optional[str]
    experiment_id: UUID
    protein_id: UUID
    ligand_id: UUID
    pocket_ids: List[int]


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool
