from __future__ import annotations

import dataclasses
from enum import Enum
from typing import List

import pydantic

from conformations.mixins import BaseModelMixin


class Integrators(Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


class GromacsForceFields(Enum):
    AMBER03 = 'amber03'
    AMBER94 = 'amber94'
    AMBER96 = 'amber96'
    AMBER99 = 'amber99'
    AMBER99SB_ILDN = 'amber99db-ildn'
    AMBER99SB = 'amber99sb'
    AMBERGS = 'amberGS'
    CHARMM27 = 'charmm27'


class GromacsWaterForceFields(Enum):
    SPC = 'spc'
    SPCE = 'spce'
    IP3P = 'ip3p'
    TIP4P = 'tip4p'
    TIP5P = 'tip5p'
    TIPS3P = 'tips3p'


class OpenMmForceFields(Enum):
    AMBER14_ALL = 'amber14-all.xml'
    AMBER14_FF14SB = 'amber14/protein.ff14SB.xml'
    AMBER14_FF15IPQ = 'amber14/protein.ff15ipq.xml'
    AMBER14_DNA_OL15 = 'amber14/DNA.OL15.xml'
    AMBER14_DNA_BSC1 = 'amber14/DNA.bsc1.xml'
    AMBER14_RNA_OL3 = 'amber14/RNA.OL3.xml'
    AMBER14_LIPID17 = 'amber14/lipid17.xml'
    AMBER14_GLYCAM_06J_1 = 'amber14/GLYCAM_06j-1.xml'
    CHARMM36 = 'charmm36.xml'


class OpenMmWaterForceFields(Enum):
    AMBER14_TIP3P = 'amber14/tip3p.xml'
    AMBER14_TIP3PFB = 'amber14/tip3pfb.xml'
    AMBER14_TIP4PEW = 'amber14/tip4pew.xml'
    AMBER14_TIP4PFB = 'amber14/tip4pfb.xml'
    AMBER14_SPCE = 'amber14/spce.xml'
    AMBER14_OPC = 'amber14/opc.xml'
    AMBER14_OPC3 = 'amber14/opc3.xml'
    CHARMM36_WATER = 'charmm36/water.xml'
    CHARMM36_SPCE = 'charmm36/spce.xml'
    CHARMM36_TIP3P_PME_B = 'charmm36/tip3p-pme-b.xml'
    CHARMM36_TIP3P_PME_F = 'charmm36/tip3p-pme-f.xml'
    CHARMM36_TIP4PEW = 'charmm36/tip4pew.xml'
    CHARMM36_TIP4P2005 = 'charmm36/tip4p2005.xml'
    CHARMM36_TIP5P = 'charmm36/tip5p.xml'
    CHARMM36_TIP5PEW = 'charmm36/tip5pew.xml'


@pydantic.dataclasses.dataclass
class GenGromacsFilesRequest(BaseModelMixin):
    pdb_content: str
    force_field: GromacsForceFields
    water_force_field: GromacsWaterForceFields


@pydantic.dataclasses.dataclass
class RunGromacsSimulationsRequest(BaseModelMixin):
    top: str
    gro: str
    temperature_k: float = 273.15
    friction_coeff: float = 1.0
    step_size: float = 0.002
    integrator: Integrators = Integrators.langevin
    take_frame_every: int = 1000
    total_frames: int = 10000


@pydantic.dataclasses.dataclass
class RunPdbSimulationsRequest(BaseModelMixin):
    pdb_content: str
    force_field: OpenMmForceFields
    water_force_field: OpenMmWaterForceFields
    temperature_k: float = 273.15
    friction_coeff: float = 1.0
    step_size: float = 0.002
    integrator: Integrators = Integrators.langevin
    take_frame_every: int = 1000
    total_frames: int = 10000


@pydantic.dataclasses.dataclass
class RunPdbFixerRequest(BaseModelMixin):
    pdb_content: str
    replace_nonstandard_residues: bool = False
    add_missing_atoms: bool = False
    add_missing_hydrogens: bool = True


@pydantic.dataclasses.dataclass
class RunPdbFixerResponse(BaseModelMixin):
    errors: List[str]
    pdb_content: str | None = None


@pydantic.dataclasses.dataclass
class RunSimulationsResponse(BaseModelMixin):
    errors: List[str]
    pdb_content: str | None = None


@pydantic.dataclasses.dataclass
class GenGroTopRequest(BaseModelMixin):
    force_field: GromacsForceFields
    water_force_field: GromacsWaterForceFields
    pdb_content: str
    ignore_missing_atoms: bool = False


@pydantic.dataclasses.dataclass
class GenGroTopResponse(BaseModelMixin):
    gro: str | None
    top: str | None
    errors: List[str]
