from pydantic import dataclasses as pcdataclass, model_validator
import datetime
from typing import List, Optional, Any

from fastapi import UploadFile, File

from nolabs.api_models.experiment import ExperimentMetadataResponse





@pcdataclass.dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    cytosolic_proteins: float
    mitochondrial_proteins: float
    nuclear_proteins: float
    other_proteins: float
    extracellular_secreted_proteins: float


@pcdataclass.dataclass
class RunLocalisationResponse:
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
