from __future__ import annotations

import dataclasses
from typing import List, Optional, Dict, Union

from chembl_query.mixins import BaseModelMixin


@dataclasses.dataclass
class Molecule:
    chembl_id: str
    molecule_type: str
    pref_name: Optional[str]
    synonyms: List[str]
    smiles: str
    link: str


@dataclasses.dataclass
class ChEMBLMoleculeRequest(BaseModelMixin):
    filters: Dict[str, Union[str, bool, int]] = None
    order_by: str = None
    limit: int = 20
    job_id: Optional[str] = None


@dataclasses.dataclass
class ChEMBLMoleculeResponse(BaseModelMixin):
    molecules: List[Molecule]


@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool


@dataclasses.dataclass
class DrugIndicationRequest(BaseModelMixin):
    filters: Dict[str, Union[str, bool, int]] = None
    order_by: str = None
    limit: int = 20  # Default page size
    job_id: Optional[str] = None


@dataclasses.dataclass
class DrugIndicationResponse(BaseModelMixin):
    drugs: List[Molecule]  # Assuming Molecule can also represent a drug in this context
    total_count: int
    page: int
    pages: int
