from pydantic import dataclasses as pcdataclass
import datetime
from typing import List

import pydantic
from fastapi import UploadFile


@pcdataclass.dataclass
class RunSolubilityRequest:
    experimentName: str
    experimentId: str
    aminoAcidSequence: str | None = None
    fasta: UploadFile | None = None


@pcdataclass.dataclass
class RunSolubilityResponse:
    solubleProbability: float | None = None
    errors: List[str] = pcdataclass.Field(default_factory=list)


@pcdataclass.dataclass
class GetExperimentRequest:
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
    solubleProbability: float


@pcdataclass.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str

@pcdataclass.dataclass
class GenerateUuidResponse:
    uuid: str
