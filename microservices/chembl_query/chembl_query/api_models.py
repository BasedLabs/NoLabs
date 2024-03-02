from __future__ import annotations

import dataclasses
from typing import List, Optional

from chembl_query.mixins import BaseModelMixin

@dataclasses.dataclass
class Molecule:
    chembl_id: str
    molecule_type: str
    pref_name: Optional[str]
    synonyms: List[str]
    smiles: str

@dataclasses.dataclass
class ChEMBLMoleculeRequest(BaseModelMixin):
    search_term: str
    job_id: Optional[str] = None
@dataclasses.dataclass
class ChEMBLMoleculeResponse(BaseModelMixin):
    molecules: List[Molecule]

@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool
