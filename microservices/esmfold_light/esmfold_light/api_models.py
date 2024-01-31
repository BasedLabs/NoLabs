from __future__ import annotations

import dataclasses


@dataclasses.dataclass
class RunEsmFoldPredictionRequest:
    protein_sequence: str
    job_id: str = None


@dataclasses.dataclass
class RunEsmFoldPredictionResponse:
    pdb_content: str
