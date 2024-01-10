from pydantic import dataclasses as pcdataclass
import datetime
from typing import List, Optional

from fastapi import UploadFile
from pydantic.dataclasses import dataclass


@pcdataclass.dataclass
class RunProteinDesignRequest:
    experimentName: str
    experimentId: str
    pdbFile: UploadFile
    contig: str = '50'
    numberOfDesign: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[int] = None


@pcdataclass.dataclass
class RunProteinDesignResponse:
    pdbContents: List[str] = pcdataclass.dataclass.field(default_factory=list)
    errors: List[str] = pcdataclass.dataclass.field(default_factory=list)


@pcdataclass.dataclass
class GetResultsRequest:
    experimentId: str


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
    metaData: ExperimentMetadataResponse
    data: str


@pcdataclass.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str


@pcdataclass.dataclass
class GenerateUuidResponse:
    uuid: str
