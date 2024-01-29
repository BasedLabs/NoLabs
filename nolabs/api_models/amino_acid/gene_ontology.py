from pydantic import dataclasses as pcdataclass, model_validator
from typing import List, Dict, Optional, Any

from fastapi import UploadFile


@pcdataclass.dataclass
class RunGeneOntologyResponseDataNode:
    name: str
    namespace: str
    edges: Dict[str, List[str]]


@pcdataclass.dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    go: Dict[str, RunGeneOntologyResponseDataNode]


@pcdataclass.dataclass
class RunGeneOntologyResponse:
    experiment_id: str
    experiment_name: str
    amino_acids: List[AminoAcidResponse]


@pcdataclass.dataclass
class ExperimentFastaPropertyResponse:
    filename: str
    content: str


@pcdataclass.dataclass
class ExperimentPropertiesResponse:
    amino_acid_sequence: str | None
    fastas: List[ExperimentFastaPropertyResponse]


@pcdataclass.dataclass
class GetExperimentResponse:
    experiment_id: str
    experiment_name: str
    amino_acids: List[AminoAcidResponse]
    properties: ExperimentPropertiesResponse
