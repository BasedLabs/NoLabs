from pydantic import dataclasses as pcdataclass
import datetime
from enum import Enum
from typing import List

from fastapi import UploadFile


class IntegratorsRequest(Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


@pcdataclass.dataclass
class RunSimulationsRequest:
    pdbFile: UploadFile
    experimentName: str
    experimentId: str
    totalFrames: int = 10000
    temperatureK: float = 273.15
    takeFrameEvery: int = 1000
    stepSize: float = 0.002
    replaceNonStandardResidues: bool = False
    addMissingAtoms: bool = False
    addMissingHydrogens: bool = True
    frictionCoeff: float = 1.0
    ignoreMissingAtoms: bool = False
    integrator: IntegratorsRequest = IntegratorsRequest.langevin


@pcdataclass.dataclass
class RunSimulationsResponse:
    pdbContent: str | None = None
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
    data: str


@pcdataclass.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str


@pcdataclass.dataclass
class GenerateUuidResponse:
    uuid: str
