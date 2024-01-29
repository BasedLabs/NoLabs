from __future__ import annotations

import dataclasses
from typing import List, Optional

import pydantic

from umol.mixins import BaseModelMixin, ErrorResponseMixin


@pydantic.dataclasses.dataclass
@dataclasses.dataclass
class RunUmolPredictionRequest(BaseModelMixin):
    protein_sequence: str
    ligand_smiles: str
    msa_content: str
    pocket_ids: List[int]
    job_id: str = None


@pydantic.dataclasses.dataclass
@dataclasses.dataclass
class RunUmolPredictionResponse(BaseModelMixin, ErrorResponseMixin):
    sdf_contents: Optional[str] = None
    pdb_contents: Optional[str] = None
    plddt_array: List[int] = dataclasses.field(default_factory=list)

@pydantic.dataclasses.dataclass
@dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool
