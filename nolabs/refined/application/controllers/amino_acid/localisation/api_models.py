from typing import List
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    cytosolic: float
    mitochondrial: float
    nuclear: float
    other: float
    extracellular: float


@dataclass
class RunJobResponse:
    job_id: UUID
    amino_acids: List[AminoAcidResponse]


@dataclass
class JobFastaPropertyResponse:
    filename: str
    content: str


@dataclass
class JobPropertiesResponse:
    fastas: List[JobFastaPropertyResponse]


@dataclass
class GetJobMetadataResponse:
    job_id: UUID
    job_name: str


@dataclass
class GetJobResponse:
    job_id: UUID
    job_name: str
    amino_acids: List[AminoAcidResponse]
    properties: JobPropertiesResponse


@dataclass
class UpdateJobRequest:
    job_name: str
