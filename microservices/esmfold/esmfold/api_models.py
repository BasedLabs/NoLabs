from __future__ import annotations

import dataclasses
from pydantic import BaseModel
from esmfold.mixins import BaseModelMixin


@dataclasses.dataclass(kw_only=True)
class RunEsmFoldPredictionRequest(BaseModelMixin, BaseModel):
    protein_sequence: str


@dataclasses.dataclass(kw_only=True)
class RunEsmFoldPredictionResponse(BaseModelMixin, BaseModel):
    pdb_content: str
