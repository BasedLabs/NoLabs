from typing import List, Optional, Any
from uuid import UUID

from pydantic import BaseModel, model_validator
from pydantic.dataclasses import dataclass


@dataclass
class JobResult:
    protein_id: UUID
    soluble_probability: float


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    proteins: List[UUID]
    result: List[JobResult]


@dataclass
class SetupJobRequest:
    experiment_id: UUID
    proteins: List[UUID]

    job_id: Optional[UUID] = None
    job_name: Optional[str] = None

    @classmethod
    @model_validator(mode='after')
    def check_inputs(cls, data: Any) -> Any:
        if not data.amino_acids:
            raise ValueError('You did not provide proteins')
        return data


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool