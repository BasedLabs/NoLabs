from __future__ import annotations

from enum import Enum
from typing import List

import pydantic
from pydantic import Field

from microservice.mixins import BaseModelMixin


@pydantic.dataclasses.dataclass
class RunRosettaFoldRequest(BaseModelMixin):
    fasta_content: str | None
    a3m_content: str | None


@pydantic.dataclasses.dataclass
class RunRosettaFoldResponse(BaseModelMixin):
    pdb_content: str | None
    errors: List[str]