from __future__ import annotations

import dataclasses
from typing import List, Optional

from esmfold.mixins import BaseModelMixin


@dataclasses.dataclass
class RunEsmFoldPredictionRequest(BaseModelMixin):
    protein_sequence: str
    job_id: str = None

@dataclasses.dataclass
class RunEsmFoldPredictionResponse(BaseModelMixin):
    errors: List[str]
    pdb_content: Optional[str] = None
