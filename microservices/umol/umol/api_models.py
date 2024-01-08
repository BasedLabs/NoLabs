from __future__ import annotations

import dataclasses
from typing import List, Optional

import pydantic
from fastapi import UploadFile

from umol.mixins import BaseModelMixin, ErrorResponseMixin


@pydantic.dataclasses.dataclass
@dataclasses.dataclass
class RunUmolPredictionRequest(BaseModelMixin):
    protein_sequence: str
    ligand_smiles: str
    msa_content: str
    pocket_ids: List[int]


@pydantic.dataclasses.dataclass
@dataclasses.dataclass
class RunUmolPredictionResponse(BaseModelMixin, ErrorResponseMixin):
    sdf_contents: Optional[str] = None
    plddt_array: List[int] = dataclasses.field(default_factory=list)
