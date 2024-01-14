from pydantic import dataclasses as pcdataclass
import datetime
from enum import Enum
from typing import List, Optional

from fastapi import UploadFile


class IntegratorsRequest(Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


@pcdataclass.dataclass
class RunSimulationsRequest:
    pdb_file: UploadFile
    experiment_name: str
    experiment_id: Optional[str]
    total_frames: Optional[int] = 10000
    temperature_k: Optional[float] = 273.15
    take_frame_every: Optional[int] = 1000
    step_size: Optional[float] = 0.002
    replace_non_standard_residues: Optional[bool] = False
    add_missing_atoms: Optional[bool] = False
    add_missing_hydrogens: Optional[bool] = True
    friction_coeff: Optional[float] = 1.0
    ignore_missing_atoms: Optional[bool] = False
    integrator: Optional[IntegratorsRequest] = IntegratorsRequest.langevin


@pcdataclass.dataclass
class RunSimulationsResponse:
    experiment_id: str
    pdb_content: str | None = None
    errors: List[str] = pcdataclass.Field(default_factory=list)


@pcdataclass.dataclass
class GetExperimentRequest:
    experiment_id: str


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
    metadata: ExperimentMetadataResponse
    data: str


@pcdataclass.dataclass
class ChangeExperimentNameRequest:
    id: str
    name: str


@pcdataclass.dataclass
class GenerateUuidResponse:
    uuid: str
