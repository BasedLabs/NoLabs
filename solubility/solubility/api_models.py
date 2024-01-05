import dataclasses
from typing import Optional

import pydantic.dataclasses

from solubility.mixins import BaseModelMixin, ErrorResponseMixing


@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunSolubilityPredictionRequest(BaseModelMixin):
    proteinSequence: str


@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunSolubilityPredictionResponse(BaseModelMixin, ErrorResponseMixing):
    solubleConfidence: Optional[float]
