from pydantic import dataclasses as pcdataclass, model_validator
import datetime
from typing import List, Optional, Any

from fastapi import UploadFile, File


@pcdataclass.dataclass
class RunSolubilityRequest:
    experiment_name: str
    experiment_id: Optional[str]
    amino_acid_sequence: Optional[str]
    fastas: Optional[List[UploadFile]]

    @model_validator(mode='after')
    @classmethod
    def check_inputs(cls, data: Any) -> Any:
        if not isinstance(data, RunSolubilityRequest):
            raise ValueError('Incorrect data type')
        if not data.amino_acid_sequence and not data.fastas:
            raise ValueError('Either specify aminoacid sequence or fastas files')
        return data


@pcdataclass.dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    soluble_probability: float


@pcdataclass.dataclass
class RunSolubilityResponse:
    experiment_id: str
    amino_acids: List[AminoAcidResponse]
    errors: List[str] = pcdataclass.Field(default_factory=list)


@pcdataclass.dataclass
class GetExperimentRequest:
    experiment_id: str


@pcdataclass.dataclass
class DeleteExperimentRequest:
    id: str


@pcdataclass.dataclass
class ExperimentMetadataResponse:
    id: str
    name: str
    date: datetime.datetime


@pcdataclass.dataclass
class GetExperimentResponse:
    metadata: ExperimentMetadataResponse
    amino_acids: List[AminoAcidResponse]


@pcdataclass.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str

@pcdataclass.dataclass
class GenerateUuidResponse:
    uuid: str
