from typing import List
from uuid import UUID

import pydantic.dataclasses
from gene_ontology.mixins import BaseModelMixin


@pydantic.dataclasses.dataclass
class RunGeneOntologyPredictionRequest(BaseModelMixin):
    job_id: UUID
    amino_acid_sequence: str


@pydantic.dataclasses.dataclass
class GoConfidenceResponse(BaseModelMixin):
    name: str
    confidence: float


@pydantic.dataclasses.dataclass
class RunGeneOntologyPredictionResponse(BaseModelMixin):
    go_confidence: List[GoConfidenceResponse]
    errors: List[str]

@pydantic.dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool