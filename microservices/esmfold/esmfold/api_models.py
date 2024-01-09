from __future__ import annotations

import dataclasses
from pydantic import BaseModel
from esmfold.mixins import BaseModelMixin


@dataclasses.dataclass
class RunEsmFoldPredictionRequest:
    protein_sequence: str

@dataclasses.dataclass
class RunEsmFoldPredictionResponse:
    pdb_content: str
