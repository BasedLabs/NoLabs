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
    experiment_id: str
    integrator: IntegratorsRequest = IntegratorsRequest.langevin
    total_frames: Optional[int] = 10000
    temperature_k: Optional[float] = 273.15
    take_frame_every: Optional[int] = 1000
    step_size: Optional[float] = 0.002
    replace_non_standard_residues: Optional[bool] = False
    add_missing_atoms: Optional[bool] = False
    add_missing_hydrogens: Optional[bool] = True
    friction_coeff: Optional[float] = 1.0
    ignore_missing_atoms: Optional[bool] = False


@pcdataclass.dataclass
class TimelineResponse:
    message: str
    error: str | None
    created_at: datetime.datetime


@pcdataclass.dataclass
class RunSimulationsResponse:
    experiment_id: str
    experiment_name: str
    timeline: List[TimelineResponse]
    pdb_content: str | None = None


@pcdataclass.dataclass
class GetExperimentRequest:
    experiment_id: str


@pcdataclass.dataclass
class DeleteExperimentRequest:
    id: str


@pcdataclass.dataclass
class ExperimentPropertiesResponse:
    pdb_file: str
    pdb_file_name: str
    total_frames: int
    temperature_k: float
    take_frame_every: int
    step_size: float
    replace_non_standard_residues: bool
    add_missing_atoms: bool
    add_missing_hydrogens: bool
    friction_coeff: float
    ignore_missing_atoms: bool
    integrator: IntegratorsRequest = IntegratorsRequest.langevin


@pcdataclass.dataclass
class GetExperimentResponse:
    experiment_id: str
    experiment_name: str
    pdb_file: str
    properties: ExperimentPropertiesResponse
