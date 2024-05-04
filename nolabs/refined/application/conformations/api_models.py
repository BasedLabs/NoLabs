from enum import Enum
from typing import List, Optional
from uuid import UUID

from fastapi import UploadFile
from pydantic.dataclasses import dataclass

from nolabs.api_models.experiment import TimelineResponse


class IntegratorsRequest(Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


@dataclass
class RunSimulationsRequest:
    pdb_file: UploadFile
    experiment_id: UUID
    job_id: Optional[UUID]
    integrator: IntegratorsRequest = IntegratorsRequest.langevin
    total_frames: int = 10000
    temperature_k: float = 273.15
    take_frame_every: int = 1000
    step_size: float = 0.002
    replace_non_standard_residues: bool = False
    add_missing_atoms: bool = False
    add_missing_hydrogens: bool = True
    friction_coeff: float = 1.0
    ignore_missing_atoms: bool = False


@dataclass
class RunSimulationsResponse:
    job_id: UUID
    timeline: List[TimelineResponse]
    pdb_content: str | None = None


@dataclass
class JobPropertiesResponse:
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


@dataclass
class GetJobResponse:
    job_id: UUID
    job_name: str
    properties: JobPropertiesResponse
    timeline: List[TimelineResponse]
    pdb_file: str | None
