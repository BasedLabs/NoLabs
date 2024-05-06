from __future__ import annotations

from pydantic.dataclasses import dataclass


@dataclass
class RunEsmFoldPredictionRequest:
    protein_sequence: str
    job_id: str = None


@dataclass
class RunEsmFoldPredictionResponse:
    pdb_content: str

@dataclass
class IsJobRunningResponse:
    is_running: bool