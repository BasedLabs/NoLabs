from typing import List
from uuid import UUID

from pydantic.dataclasses import dataclass

from localisation.mixins import BaseModelMixin


@dataclass
class RunLocalisationPredictionRequest(BaseModelMixin):
    job_id: UUID
    amino_acid_sequence: str


@dataclass
class RunLocalisationPredictionResponse(BaseModelMixin):
    cytosolic_proteins: float
    mitochondrial_proteins: float
    nuclear_proteins: float
    other_proteins: float
    extracellular_secreted_proteins: float
    errors: List[str]

@dataclass
class IsJobRunningResponse:
    is_running: bool