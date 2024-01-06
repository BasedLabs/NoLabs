from __future__ import annotations

import dataclasses
from typing import List

from pydantic import BaseModel, Field

from protein_design.mixins import BaseModelMixin, ErrorResponseMixing


@dataclasses.dataclass(kw_only=True)
class RunRfdiffusionRequest(BaseModelMixin, BaseModel):
    pdbContent: str
    contig: str = '5'
    timesteps: int = 10
    hotspots: str | None
    numberOfDesigns: int = 1


@dataclasses.dataclass(kw_only=True)
class RunRfdiffusionResponse(BaseModelMixin, ErrorResponseMixing, BaseModel):
    pdbsContents: List[str] = Field(default_factory=list)
