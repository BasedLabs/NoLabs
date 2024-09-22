from typing import Any, List, Optional
from uuid import UUID

from pydantic import model_validator
from pydantic.dataclasses import dataclass


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    protein_id: UUID
    result: List[int]


@dataclass
class SetupJobRequest:
    experiment_id: UUID

    protein_id: UUID

    job_id: Optional[UUID] = None
    job_name: Optional[str] = None

    @classmethod
    @model_validator(mode="after")
    def check_inputs(cls, data: Any) -> Any:
        if not data.amino_acids:
            raise ValueError("You did not provide proteins")
        return data


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool
    result_valid: bool
