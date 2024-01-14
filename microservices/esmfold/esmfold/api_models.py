from __future__ import annotations

import dataclasses
from esmfold.mixins import BaseModelMixin


@dataclasses.dataclass
class RunEsmFoldPredictionRequest(BaseModelMixin):
    protein_sequence: str

@dataclasses.dataclass
class RunEsmFoldPredictionResponse(BaseModelMixin):
    pdb_content: str
