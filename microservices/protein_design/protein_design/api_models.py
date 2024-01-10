from __future__ import annotations

import dataclasses
from typing import List

import pydantic.dataclasses
from pydantic import Field

from protein_design.mixins import BaseModelMixin, ErrorResponseMixing


@pydantic.dataclasses.dataclass
@dataclasses.dataclass()
class RunRfdiffusionRequest(BaseModelMixin):
    pdbContent: str
    hotspots: str | None
    contig: str = '5'
    timesteps: int = 10
    numberOfDesigns: int = 1


@pydantic.dataclasses.dataclass
@dataclasses.dataclass()
class RunRfdiffusionResponse(BaseModelMixin, ErrorResponseMixing):
    pdbsContents: List[str] = Field(default_factory=list)
