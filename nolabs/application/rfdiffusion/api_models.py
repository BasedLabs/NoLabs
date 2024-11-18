from __future__ import annotations
from typing import List, Optional
from uuid import UUID

from pydantic import BaseModel


class JobResponse(BaseModel):
    job_id: UUID
    job_name: str
    protein_id: UUID
    binder_ids: List[UUID]
    contig: Optional[str] = None
    number_of_designs: Optional[int] = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None


class SetupJobRequest(BaseModel):
    protein_id: UUID
    contig: str
    number_of_designs: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None
    job_id: Optional[UUID] = None
    job_name: Optional[str] = None
