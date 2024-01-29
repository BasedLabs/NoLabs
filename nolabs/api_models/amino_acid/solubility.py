from pydantic import dataclasses as pcdataclass, model_validator
import datetime
from typing import List, Optional, Any

from fastapi import UploadFile, File

from nolabs.api_models.experiment import ExperimentMetadataResponse



@pcdataclass.dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    soluble_probability: float


@pcdataclass.dataclass
class RunSolubilityResponse:
    experiment_id: str
    experiment_name: str
    amino_acids: List[AminoAcidResponse]
    errors: List[str] = pcdataclass.Field(default_factory=list)


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
