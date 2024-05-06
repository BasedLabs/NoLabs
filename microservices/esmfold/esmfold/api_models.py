from __future__ import annotations

from pydantic.dataclasses import dataclass
from typing import List, Optional

from esmfold.mixins import BaseModelMixin


@dataclass
class RunEsmFoldPredictionRequest(BaseModelMixin):
    protein_sequence: str
    job_id: str = None

@dataclass
class RunEsmFoldPredictionResponse(BaseModelMixin):
    errors: List[str]
    pdb_content: Optional[str] = None

@dataclass
class IsJobRunningResponse:
    is_running: bool
