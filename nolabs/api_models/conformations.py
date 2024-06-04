from pydantic import dataclasses as pcdataclass
import datetime
from enum import Enum
from typing import List, Optional

from fastapi import UploadFile

from nolabs.api_models.experiment import TimelineResponse


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
    total_frames: int = 10000
    temperature_k: float = 309.75
    take_frame_every: int = 1000
    step_size: float = 0.002
    replace_non_standard_residues: bool = False
    add_missing_atoms: bool = False
    add_missing_hydrogens: bool = True
    friction_coeff: float = 1.0
    ignore_missing_atoms: bool = False


@pcdataclass.dataclass
class RunSimulationsResponse:
    experiment_id: str
    experiment_name: str
    timeline: List[TimelineResponse]
    pdb_content: str | None = None


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
    properties: ExperimentPropertiesResponse
    timeline: List[TimelineResponse]
    pdb_file: str | None
