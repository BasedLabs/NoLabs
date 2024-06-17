from enum import Enum
from typing import List, Optional
from uuid import UUID

from pydantic.dataclasses import dataclass

from nolabs.refined.application.use_cases.experiments.api_models import TimelineResponse


class IntegratorsRequest(Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


@dataclass
class JobResponse:
    job_id: UUID
    job_name: str
    experiment_id: UUID
    protein_id: UUID

    timeline: List[TimelineResponse]

    total_frames: int
    temperature_k: float
    take_frame_every: int
    step_size: float
    replace_non_standard_residues: bool
    add_missing_atoms: bool
    add_missing_hydrogens: bool
    friction_coeff: float
    ignore_missing_atoms: bool
    integrator: IntegratorsRequest

    md_content: Optional[bytes] = None


@dataclass
class SetupJobRequest:
    protein_id: UUID

    job_id: Optional[UUID] = None
    job_name: Optional[str] = None

    total_frames: int = 10000
    temperature_k: float = 273.15
    take_frame_every: int = 1000
    step_size: float = 0.002
    replace_non_standard_residues: bool = False
    add_missing_atoms: bool = False
    add_missing_hydrogens: bool = True
    friction_coeff: float = 1.0
    ignore_missing_atoms: bool = False
    integrator: IntegratorsRequest = IntegratorsRequest.langevin

    experiment_id: Optional[UUID] = None


@dataclass
class RunJobRequest:
    job_id: UUID


@dataclass
class GetJobStatusResponse:
    running: bool
