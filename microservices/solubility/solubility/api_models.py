import dataclasses
from typing import Optional

import pydantic.dataclasses

from solubility.mixins import BaseModelMixin, ErrorResponseMixing


@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunSolubilityPredictionRequest(BaseModelMixin):
    amino_acid_sequence: str


@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunSolubilityPredictionResponse(BaseModelMixin, ErrorResponseMixing):
    soluble_probability: Optional[float]
