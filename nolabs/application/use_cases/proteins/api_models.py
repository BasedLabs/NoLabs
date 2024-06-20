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
class ProteinMetadataResponse:
    id: UUID
    name: str
    experiment_id: UUID

    binding_pockets: List[int]
    fasta_name: str
    pdb_name: str
    localisation: Optional[ProteinLocalisationResponse] = None
    gene_ontology: Optional[Dict[str, Any]] = None
    soluble_probability: Optional[float] = None
    link: Optional[str] = None


@dataclass
class ProteinContentResponse:
    id: UUID
    name: str
    experiment_id: UUID

    binding_pockets: List[int]
    fasta_name: str
    pdb_name: str
    localisation: Optional[ProteinLocalisationResponse] = None
    gene_ontology: Optional[Dict[str, Any]] = None
    soluble_probability: Optional[float] = None
    msa: Optional[str] = None
    md_pdb_content: Optional[str] = None
    fasta_content: Optional[str] = None
    pdb_content: Optional[str] = None
    link: Optional[str] = None

@dataclass
class ProteinSearchQuery:
    name: Optional[str] = None
    experiment_id: Optional[UUID] = None
    ids: Optional[List[UUID]] = None


@dataclass
class ProteinSearchMetadataQuery:
    name: Optional[str] = None
    experiment_id: Optional[UUID] = None
    ids: Optional[List[UUID]] = None



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