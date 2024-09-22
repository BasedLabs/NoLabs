from __future__ import annotations

from typing import List, Optional
from uuid import UUID

import pydantic
from protein_design.mixins import BaseModelMixin
from pydantic import Field
from pydantic.dataclasses import dataclass


@dataclass
class RunRfdiffusionRequest(BaseModelMixin):
    job_id: UUID
    pdb_content: str
    contig: str
    hotspots: Optional[str] = ""
    timesteps: Optional[int] = 10
    number_of_designs: Optional[int] = 1


@dataclass
class RunRfdiffusionResponse(BaseModelMixin):
    pdbs_content: List[str] = Field(default_factory=list)
    errors: List[str] = pydantic.dataclasses.Field(default_factory=list)


@dataclass
class IsJobRunningResponse:
    is_running: bool
