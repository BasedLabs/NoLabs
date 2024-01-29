from pydantic import dataclasses as pcdataclass, model_validator
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
