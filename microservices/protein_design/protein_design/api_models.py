from __future__ import annotations

import dataclasses
from typing import List, Optional

import pydantic.dataclasses
from pydantic import Field

from protein_design.mixins import BaseModelMixin, ErrorResponseMixing


@pydantic.dataclasses.dataclass
@dataclasses.dataclass()
class RunRfdiffusionRequest(BaseModelMixin):
    pdb_content: str
    hotspots: Optional[str] = ''
    contig: Optional[str] = '5'
    timesteps: Optional[int] = 10
    number_of_designs: Optional[int] = 1


@pydantic.dataclasses.dataclass
@dataclasses.dataclass()
class RunRfdiffusionResponse(BaseModelMixin, ErrorResponseMixing):
    pdbs_content: List[str] = Field(default_factory=list)
