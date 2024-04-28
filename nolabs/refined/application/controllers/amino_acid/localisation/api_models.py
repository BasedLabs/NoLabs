from typing import List
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    cytosolic_proteins: float
    mitochondrial_proteins: float
    nuclear_proteins: float
    other_proteins: float
    extracellular_secreted_proteins: float


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
class GetJobResponse:
    job_id: UUID
    job_name: str
    amino_acids: List[AminoAcidResponse]
    properties: JobPropertiesResponse