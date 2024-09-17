from __future__ import annotations

__all__ = [
    "RunReinforcementLearningRequest",
    "RunSamplingRequest",
    "PreparePdbqtRequest",
    "PreparePdbqtResponse",
]

from pydantic import BaseModel


class RunReinforcementLearningRequest(BaseModel):
    config_id: str


class RunSamplingRequest(BaseModel):
    config_id: str
    number_of_molecules_to_generate: int


class PreparePdbqtRequest(BaseModel):
    pdb: bytes


class PreparePdbqtResponse(BaseModel):
    pdbqt: bytes
