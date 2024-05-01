from typing import List, Dict
from uuid import UUID

from pydantic.dataclasses import dataclass


@dataclass
class RunGeneOntologyResponseDataNode:
    name: str
    namespace: str
    edges: Dict[str, List[str]]


@dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    go: Dict[str, RunGeneOntologyResponseDataNode]


@dataclass
class RunGeneOntologyResponse:
    experiment_id: str
    experiment_name: str
    amino_acids: List[AminoAcidResponse]


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
