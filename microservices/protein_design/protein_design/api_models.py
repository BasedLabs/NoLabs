from __future__ import annotations

import dataclasses
from typing import List

import pydantic.dataclasses
from pydantic import Field

from protein_design.mixins import BaseModelMixin, ErrorResponseMixing


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunRfdiffusionRequest(BaseModelMixin):
    pdbContent: str
    contig: str = '5'
    timesteps: int = 10
    hotspots: str | None
    numberOfDesigns: int = 1


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunRfdiffusionResponse(BaseModelMixin, ErrorResponseMixing):
    pdbsContents: List[str] = Field(default_factory=list)
