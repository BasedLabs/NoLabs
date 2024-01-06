from __future__ import annotations

import dataclasses
from typing import List
import pydantic
from fastapi import UploadFile

from p2rank.mixins import BaseModelMixin


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunP2RankPredictionRequest(BaseModelMixin):
    protein_file: UploadFile


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunP2RankPredictionResponse(BaseModelMixin):
    pocket_ids: List[int]
