from pydantic import dataclasses as pcdataclass
import datetime
from typing import List, Dict

import pydantic
from fastapi import UploadFile


@pcdataclass.dataclass
class RunGeneOntologyRequest:
    experimentName: str
    experimentId: str
    aminoAcidSequence: str | None = None
    fasta: UploadFile | None = None


@pcdataclass.dataclass
class RunGeneOntologyResponseDataNode:
    name: str
    namespace: str
    edges: Dict[str, List[str]]

@pcdataclass.dataclass
class RunGeneOntologyResponse:
    data: Dict[str, RunGeneOntologyResponseDataNode] | None = None
    errors: List[str] = pcdataclass.dataclass.field(default_factory=list)


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
    data: Dict[str, RunGeneOntologyResponseDataNode]


@pcdataclass.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str

@pcdataclass.dataclass
class GenerateUuidResponse:
    uuid: str
