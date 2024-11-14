__all__ = ["JobResult", "JobResponse", "SetupJobRequest", "RunJobRequest"]

from typing import List, Optional, Dict
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class JobResult:
    sequence_id: UUID
    sequence: str
    fasta_content: str
    score: Optional[float] = None
    global_score: Optional[float] = None
    T: Optional[float] = None
    sample: Optional[int] = None
    seq_recovery: Optional[float] = None
    # Add additional fields if needed


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    num_seq_per_target: int
    sampling_temp: float
    seed: int
    batch_size: int
    is_homomer: bool
    protein_id: UUID
    result: List[JobResult]
    chains_to_design: Optional[List[str]] = None
    fixed_positions: Optional[Dict[str, List[int]]] = None


@dataclass
class SetupJobRequest:
    experiment_id: UUID
    protein_id: UUID
    num_seq_per_target: Optional[int] = 2
    sampling_temp: Optional[float] = 0.1
    seed: Optional[int] = 37
    batch_size: Optional[int] = 1
    is_homomer: Optional[bool] = False
    chains_to_design: Optional[List[str]] = None
    fixed_positions: Optional[Dict[str, List[int]]] = None
    job_id: Optional[UUID] = None
    job_name: Optional[str] = None


@dataclass
class RunJobRequest:
    job_id: UUID
