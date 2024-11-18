__all__ = ["JobResult", "JobResponse"]


from typing import List, Optional
from uuid import UUID

from pydantic import BaseModel


class JobResult(BaseModel):
    complex_id: UUID
    sdf_content: str
    minimized_affinity: Optional[float] = None
    scored_affinity: Optional[float] = None
    confidence: Optional[float] = None


class JobResponse(BaseModel):
    job_id: UUID
    job_name: str
    samples_per_complex: int
    protein_id: UUID
    ligand_id: UUID
    result: List[JobResult]


class SetupJobRequest(BaseModel):
    experiment_id: UUID
    protein_id: UUID
    ligand_id: UUID
    samples_per_complex: Optional[int] = 2
    job_id: Optional[UUID] = None
    job_name: Optional[str] = None
