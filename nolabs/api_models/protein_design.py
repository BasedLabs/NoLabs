import datetime
from typing import List, Optional

from fastapi import UploadFile
from pydantic import dataclasses


@dataclasses.dataclass
class RunProteinDesignRequest:
    experiment_name: str
    experiment_id: str
    pdb_file: UploadFile
    contig: str = '50'
    number_of_desings: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None


@dataclasses.dataclass
class RunProteinDesignResponse:
    experiment_id: str
    experiment_name: str
    pdbs_content: List[str] = dataclasses.Field(default_factory=list)


@dataclasses.dataclass
class GetResultsRequest:
    experimentId: str


@dataclasses.dataclass
class DeleteExperimentRequest:
    id: str


@dataclasses.dataclass
class ExperimentMetadataResponse:
    id: str
    name: str
    date: datetime.datetime


@dataclasses.dataclass
class GetExperimentResponse:
    experiment_id: str
    experiment_name: str
    pdbs_content: List[str]
    pdb_file: str
    pdb_file_name: str
    contig: str
    number_of_desings: int
    timesteps: int
    hotspots: str


@dataclasses.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str


@dataclasses.dataclass
class GenerateUuidResponse:
    uuid: str
