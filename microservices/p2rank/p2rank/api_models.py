from __future__ import annotations

import dataclasses
from typing import List
import pydantic

from p2rank.mixins import BaseModelMixin


@pydantic.dataclasses.dataclass
@dataclasses.dataclass
class RunP2RankPredictionRequest(BaseModelMixin):
    pdb_contents: str
    job_id: str = None

@pydantic.dataclasses.dataclass
@dataclasses.dataclass
class RunP2RankPredictionResponse(BaseModelMixin):
    pocket_ids: List[int]
