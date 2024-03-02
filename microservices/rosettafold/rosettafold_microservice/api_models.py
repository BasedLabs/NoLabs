from __future__ import annotations

from enum import Enum
from typing import List

import pydantic
from pydantic import Field

from rosettafold_microservice.mixins import BaseModelMixin


@pydantic.dataclasses.dataclass
class RunRosettaFoldRequest(BaseModelMixin):
    fasta_content: str = Field(min_length=1)


@pydantic.dataclasses.dataclass
class RunRosettaFoldResponse(BaseModelMixin):
    pdb_content: str | None
    errors: List[str]