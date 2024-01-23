from typing import Optional

import pydantic

from nolabs.api_models.conformations import IntegratorsRequest


@pydantic.dataclasses.dataclass
class ExperimentProperties:
    input_file_name: str
    input_file: str
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