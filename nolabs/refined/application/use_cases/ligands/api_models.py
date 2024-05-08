__all__ = [
    'LigandSearchQuery',
    'LigandResponse',
    'UploadLigandRequest'
]

from typing import Dict, Any, Optional
from uuid import UUID

from fastapi import UploadFile
from pydantic.dataclasses import dataclass


@dataclass
class LigandResponse:
    id: UUID
    name: str
    experiment_id: UUID
    smiles_content: Optional[str]
    sdf_content: Optional[str]
    drug_likeness: Optional[float]
    designed_ligand_score: Optional[float]


@dataclass
class LigandSearchQuery:
    name: Optional[str]
    experiment_id: Optional[UUID]


@dataclass
class UploadLigandRequest:
    experiment_id: UUID
    name: Optional[str]
    smiles: Optional[UploadFile]
    sdf: Optional[UploadFile]