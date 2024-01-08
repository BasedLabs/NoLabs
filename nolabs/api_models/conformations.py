import dataclasses
import datetime
from enum import Enum
from typing import List

from fastapi import UploadFile
from pydantic.dataclasses import dataclass


class IntegratorsRequest(Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


@dataclass
@dataclasses.dataclass
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


@dataclass
@dataclasses.dataclass
class RunSimulationsResponse:
    pdbContent: str | None = None
    errors: List[str] = dataclasses.field(default_factory=list)


@dataclass
@dataclasses.dataclass
class GetExperimentRequest:
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
