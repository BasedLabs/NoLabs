from typing import Dict, Any
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class ProteinLocalisation:
    cytosolic: float
    mitochondrial: float
    nuclear: float
    other: float
    extracellular: float


@dataclass
class Protein:
    id: UUID
    name: str
    experiment_id: UUID
    fasta_content: str | None
    pdb_content: str | None
    localisation: ProteinLocalisation
    gene_ontology: Dict[str, Any]
    soluble_probability: float
    msa: str | None
    md_pdb_content: str | None