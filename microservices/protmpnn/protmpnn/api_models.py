from __future__ import annotations

import dataclasses
from typing import List, Optional, Dict

from protmpnn.mixins import BaseModelMixin

@dataclasses.dataclass
class RunProtMPNNPredictionRequest(BaseModelMixin):
    pdb_contents: str  # Protein contents in PDB format
    is_homomer: bool  # Whether it's a homomer or monomer
    chains_to_design: Optional[List[str]] = None  # Chains to design (only for homomers)
    fixed_positions: Optional[Dict[str, List[int]]] = None  # Fixed positions per chain
    num_seq_per_target: int = 2
    sampling_temp: float = 0.1
    seed: int = 37
    batch_size: int = 1
    job_id: Optional[str] = None  # Optional job ID if you want to track the job status

@dataclasses.dataclass
class RunProtMPNNPredictionResponse(BaseModelMixin):
    success: bool
    message: str
    sequences: List[str]  # List of designed sequences
    fasta_contents: List[str]  # List of FASTA file contents

@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool
