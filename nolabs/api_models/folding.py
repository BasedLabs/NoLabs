import datetime
from typing import List

from pydantic import dataclasses as pcdataclass
from fastapi import UploadFile


@pcdataclass.dataclass
class RunFoldingRequest:
    experimentName: str
    experimentId: str | None
    aminoAcidSequence: str | None = None
    fasta: UploadFile | None = None


@pcdataclass.dataclass
class RunFoldingResponse:
    experimentId: str
    pdbContent: str | None = None
    errors: List[str] = pcdataclass.dataclass.field(default_factory=list)


@pcdataclass.dataclass
class GetExperimentRequest:
    experimentId: str


@pcdataclass.dataclass
class DeleteExperimentRequest:
    id: str


@pcdataclass.dataclass
class GetExperimentResponse:
    metaData: ExperimentMetadataResponse
    solubleProbability: float
