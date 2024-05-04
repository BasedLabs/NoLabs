from typing import List
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    pdb_file_name: str
    pdb_file: str


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