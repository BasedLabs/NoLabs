from typing import Optional
from uuid import UUID

from pydantic.dataclasses import dataclass

from solubility.mixins import BaseModelMixin, ErrorResponseMixing


@dataclass
class RunSolubilityPredictionRequest(BaseModelMixin):
    job_id: UUID
    amino_acid_sequence: str


@dataclass
class RunSolubilityPredictionResponse(BaseModelMixin, ErrorResponseMixing):
    soluble_probability: Optional[float]


@dataclass
class IsJobRunningResponse:
    is_running: bool
