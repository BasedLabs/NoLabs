import dataclasses
from typing import Dict, List

import pydantic.dataclasses
from gene_ontology.mixins import BaseModelMixin


@pydantic.dataclasses.dataclass
class RunGeneOntologyPredictionRequest(BaseModelMixin):
    amino_acid_sequence: str


@pydantic.dataclasses.dataclass
class GoConfidenceResponse(BaseModelMixin):
    name: str
    confidence: float


@pydantic.dataclasses.dataclass
class RunGeneOntologyPredictionResponse(BaseModelMixin):
    go_confidence: List[GoConfidenceResponse]
    errors: List[str]
