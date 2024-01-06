from __future__ import annotations

import dataclasses
from typing import List
from pydantic import BaseModel
from fastapi import UploadFile

from umol.mixins import BaseModelMixin

@dataclasses.dataclass(kw_only=True)
class RunUmolPredictionRequest(BaseModelMixin, BaseModel):
    protein_file: UploadFile
    ligand_smiles: str
    msa_file: UploadFile
    pocket_ids: List[int]

@dataclasses.dataclass(kw_only=True)
class RunUmolPredictionResponce(BaseModelMixin, BaseModel):
    sdf_contents: str
    plddt_array: List[int]