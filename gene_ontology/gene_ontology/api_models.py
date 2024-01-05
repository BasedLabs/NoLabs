import dataclasses
from typing import Dict

import pydantic.dataclasses

from gene_ontology.mixins import BaseModelMixin, ErrorResponseMixing


@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunGeneOntologyPredictionRequest(BaseModelMixin):
    proteinSequence: str


@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunGeneOntologyPredictionResponse(BaseModelMixin, ErrorResponseMixing):
    go_confidence: Dict[str, float]
