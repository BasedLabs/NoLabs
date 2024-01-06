from __future__ import annotations

import dataclasses
from typing import List
from pydantic import BaseModel
from fastapi import UploadFile

from p2rank.mixins import BaseModelMixin

@dataclasses.dataclass(kw_only=True)
class RunP2RankPredictionRequest(BaseModelMixin, BaseModel):
    protein_file: UploadFile

@dataclasses.dataclass(kw_only=True)
class RunP2RankPredictionResponse(BaseModelMixin, BaseModel):
    pocket_ids: List[int]