import dataclasses
from typing import Optional, List

import pydantic.dataclasses

from localisation.mixins import BaseModelMixin


@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunLocalisationPredictionRequest(BaseModelMixin):
    amino_acid_sequence: str


@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunLocalisationPredictionResponse(BaseModelMixin):
    cytosolic_proteins: float
    mitochondrial_proteins: float
    nuclear_proteins: float
    other_proteins: float
    extracellular_secreted_proteins: float
    errors: List[str]
