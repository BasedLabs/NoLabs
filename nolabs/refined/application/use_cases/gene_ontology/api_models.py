from typing import List, Dict, Optional, Any
from uuid import UUID

from pydantic import BaseModel, model_validator
from pydantic.dataclasses import dataclass


@dataclass
class RunGeneOntologyResponseDataNode:
    name: str
    namespace: str
    edges: Dict[str, List[str]]


@dataclass
class JobResult:
    protein_id: UUID
    go: Dict[str, RunGeneOntologyResponseDataNode]


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


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool


