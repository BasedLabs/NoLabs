from __future__ import annotations

import dataclasses
from typing import List, Optional

from diffdock.mixins import BaseModelMixin


@dataclasses.dataclass
class RunDiffDockPredictionRequest(BaseModelMixin):
    pdb_contents: str  # Protein contents in PDB format
    sdf_contents: str  # Ligand contents in SDF format
    inference_steps: int = 20
    samples_per_complex: int = 40
    batch_size: int = 10
    actual_steps: int = 18
    no_final_step_noise: bool = True
    job_id: str = None  # Optional job ID if you want to track the job status

@dataclasses.dataclass
class RunDiffDockPredictionResponse(BaseModelMixin):
    errors: List[str]
    pdb_content: Optional[str] = None

@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool
