from __future__ import annotations

import dataclasses
from typing import List, Optional

from esmfold.mixins import BaseModelMixin


@dataclasses.dataclass
class RunEsmFoldPredictionRequest(BaseModelMixin):
    protein_sequence: str

@dataclasses.dataclass
class RunEsmFoldPredictionResponse(BaseModelMixin):
    errors: List[str]
    pdb_content: Optional[str] = None
