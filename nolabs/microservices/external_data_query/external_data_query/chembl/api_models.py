from __future__ import annotations

from typing import Dict, List, Optional, Union

from external_data_query.mixins import BaseModelMixin
from pydantic.dataclasses import dataclass


@dataclass
class Molecule:
    chembl_id: str
    synonyms: List[str]
    smiles: str
    link: str
    molecule_type: str | None = None
    pref_name: str | None = None


@dataclass
class ChEMBLMoleculeRequest(BaseModelMixin):
    filters: Dict[str, Union[str, bool, int]] = None
    order_by: str = None
    limit: int = 20
    job_id: Optional[str] = None


@dataclass
class ChEMBLMoleculeResponse(BaseModelMixin):
    molecules: List[Molecule]


@dataclass
class IsJobRunningResponse:
    is_running: bool


@dataclass
class DrugIndicationRequest(BaseModelMixin):
    filters: Dict[str, Union[str, bool, int]] = None
    order_by: str = None
    limit: int = 20  # Default page size
    job_id: Optional[str] = None


@dataclass
class DrugIndicationResponse(BaseModelMixin):
    drugs: List[Molecule]  # Assuming Molecule can also represent a drug in this context
    total_count: int
    page: int
    pages: int
