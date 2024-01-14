from __future__ import annotations

import dataclasses


@dataclasses.dataclass
class RunEsmFoldPredictionRequest:
    protein_sequence: str


@dataclasses.dataclass
class RunEsmFoldPredictionResponse:
    pdb_content: str
