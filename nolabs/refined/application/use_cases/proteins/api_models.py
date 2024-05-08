from typing import Dict, Any, Optional
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
    fasta_content: Optional[str]
    pdb_content: Optional[str]
    localisation: Optional[ProteinLocalisationResponse]
    gene_ontology: Optional[Dict[str, Any]]
    soluble_probability: Optional[float]
    msa: Optional[str]
    md_pdb_content: Optional[str]


@dataclass
class ProteinSearchQuery:
    name: Optional[str]
    experiment_id: Optional[UUID]


@dataclass
class UploadProteinRequest:
    experiment_id: UUID
    name: Optional[str]
    fasta: Optional[UploadFile]
    pdb: Optional[UploadFile]