from typing import Dict, Any, Optional, List
from uuid import UUID

from fastapi import UploadFile
from pydantic.dataclasses import dataclass


@dataclass
class ProteinLocalisationResponse:
    cytosolic: float
    mitochondrial: float
    nuclear: float
    other: float
    extracellular: float


@dataclass
class ProteinResponse:
    id: UUID
    name: str
    experiment_id: UUID

    binding_pockets: List[int]
    localisation: Optional[ProteinLocalisationResponse] = None
    gene_ontology: Optional[Dict[str, Any]] = None
    soluble_probability: Optional[float] = None
    msa: Optional[str] = None
    md_pdb_content: Optional[str] = None
    fasta_content: Optional[str] = None
    pdb_content: Optional[str] = None


@dataclass
class ProteinSearchQuery:
    name: Optional[str] = None
    experiment_id: Optional[UUID] = None


@dataclass
class UploadProteinRequest:
    experiment_id: UUID
    name: Optional[str] = None
    fasta: Optional[UploadFile] = None
    pdb: Optional[UploadFile] = None


@dataclass
class UpdateProteinRequest:
    protein_id: UUID
    name: Optional[str] = None
    fasta: Optional[UploadFile] = None
    pdb: Optional[UploadFile] = None