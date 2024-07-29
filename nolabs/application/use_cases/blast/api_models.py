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
class HspModel:
    num: int
    bit_score: float
    score: int
    evalue: float
    query_from: int
    query_to: int
    hit_from: int
    hit_to: int
    query_frame: int
    hit_frame: int
    identity: int
    positive: int
    gaps: int
    align_len: int
    qseq: str
    hseq: str
    midline: str


@dataclass
class HitModel:
    num: int
    id: str
    definition: str
    accession: str
    length: int
    hsps: List[HspModel]


@dataclass
class JobResult:
    protein_id: UUID
    program: str
    database: str
    query_id: str
    query_def: str
    query_len: int
    hits: List[HitModel]


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    protein_id: UUID
    result: List[JobResult]


@dataclass
class SetupJobRequest:
    experiment_id: UUID
    protein_id: UUID
    descriptions: int = 10,
    alignments: int = 10,
    hitlist_size: int = 10,
    expect: float = 10.0,
    job_id: Optional[UUID] = None
    job_name: Optional[str] = None


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool
    result_valid: bool
