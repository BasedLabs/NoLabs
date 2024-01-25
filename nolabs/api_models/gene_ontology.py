from pydantic import dataclasses as pcdataclass, model_validator
import datetime
from typing import List, Dict, Optional, Any

from fastapi import UploadFile


@pcdataclass.dataclass
class RunGeneOntologyRequest:
    experiment_name: str
    experiment_id: Optional[str]
    amino_acid_sequence: Optional[str]
    fastas: Optional[List[UploadFile]]

    @model_validator(mode='after')
    @classmethod
    def check_inputs(cls, data: Any) -> Any:
        if not isinstance(data, RunGeneOntologyRequest):
            raise ValueError('Incorrect data type')
        if not data.amino_acid_sequence and not data.fastas:
            raise ValueError('Either specify aminoacid sequence or fastas files')
        return data


@pcdataclass.dataclass
class RunGeneOntologyResponseDataNode:
    name: str
    namespace: str
    edges: Dict[str, List[str]]


@pcdataclass.dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    go: List[RunGeneOntologyResponseDataNode]


@pcdataclass.dataclass
class RunGeneOntologyResponse:
    experiment_id: str
    amino_acids: List[AminoAcidResponse]
    errors: List[str] = pcdataclass.Field(default_factory=list)


@pcdataclass.dataclass
class GetExperimentRequest:
    experimentId: str


@pcdataclass.dataclass
class GetExperimentResponse:
    amino_acids: List[AminoAcidResponse]
