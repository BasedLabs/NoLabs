import dataclasses
import datetime
from enum import Enum
from typing import List, Optional

from fastapi import UploadFile
from pydantic.dataclasses import dataclass


@dataclass
@dataclasses.dataclass
class RunProteinDesignRequest:
    experimentName: str
    experimentId: str
    pdbFile: UploadFile
    contig: str = '50'
    numberOfDesign: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[int] = None


@dataclass
@dataclasses.dataclass
class RunProteinDesignResponse:
    pdbContents: List[str] = []
    errors: List[str] = dataclasses.field(default_factory=list)


@dataclass
@dataclasses.dataclass
class GetResultsRequest:
    experimentId: str


@dataclass
@dataclasses.dataclass
class DeleteExperimentRequest:
    id: str


@dataclass
@dataclasses.dataclass
class ExperimentMetadataResponse:
    id: str
    name: str
    date: datetime.datetime


@dataclass
@dataclasses.dataclass
class GetExperimentResponse:
    metaData: ExperimentMetadataResponse
    data: str


@dataclass
@dataclasses.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str


@dataclass
@dataclasses.dataclass
class GenerateUuidResponse:
    uuid: str
