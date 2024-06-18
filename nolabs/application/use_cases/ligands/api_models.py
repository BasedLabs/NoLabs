__all__ = [
    'LigandSearchQuery',
    'LigandResponse',
    'UploadLigandRequest'
]

from dataclasses import field
from typing import Dict, Any, Optional
from uuid import UUID

from fastapi import UploadFile
from pydantic.dataclasses import dataclass


@dataclass
class LigandResponse:
    id: UUID
    name: str
    experiment_id: UUID
    smiles_content: Optional[str] = None
    sdf_content: Optional[str] = None
    link: Optional[str] = None
    image: Optional[str] = field(default=None, repr=False)
    drug_likeness: Optional[float] = None
    designed_ligand_score: Optional[float] = None


@dataclass
class LigandSearchQuery:
    name: Optional[str] = None
    experiment_id: Optional[UUID] = None


@dataclass
class UploadLigandRequest:
    experiment_id: UUID
    name: Optional[str] = None
    smiles: Optional[UploadFile] = None
    sdf: Optional[UploadFile] = None
    link: Optional[str] = None


@dataclass
class UpdateLigandRequest:
    ligand_id: UUID
    name: Optional[str] = None
    smiles: Optional[UploadFile] = None
    sdf: Optional[UploadFile] = None
    link: Optional[str] = None