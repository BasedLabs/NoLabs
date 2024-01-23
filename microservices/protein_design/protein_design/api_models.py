from __future__ import annotations

from typing import List, Optional

import dataclasses
import pydantic.dataclasses
from pydantic import Field

from protein_design.mixins import BaseModelMixin


@pydantic.dataclasses.dataclass
class RunRfdiffusionRequest(BaseModelMixin):
    pdb_content: str
    contig: str
    hotspots: Optional[str] = ''
    timesteps: Optional[int] = 10
    number_of_designs: Optional[int] = 1


@pydantic.dataclasses.dataclass
class RunRfdiffusionResponse(BaseModelMixin):
    pdbs_content: List[str] = Field(default_factory=list)
    errors: List[str] = pydantic.dataclasses.Field(default_factory=list)
