from __future__ import annotations

import dataclasses
from enum import Enum

import pydantic
from pydantic import BaseModel

from conformations.mixins import BaseModelMixin, ErrorResponseMixing


class Integrators(Enum):
    langevin = 'LangevinIntegator'
    langevin_middle = 'LangevinMiddleIntegator'
    nose_hoover = 'NoseHooverIntegrator'
    brownian = 'BrownianIntegrator'
    variable_verlet = 'VariableVerletIntegrator'


class ForceFields(Enum):
    AMBER14_ALL = 'amber14-all.xml'
    AMBER14_FF14SB = 'amber14/protein.ff14SB.xml'
    AMBER14_FF15IPQ = 'amber14/protein.ff15ipq.xml'
    AMBER14_DNA_OL15 = 'amber14/DNA.OL15.xml'
    AMBER14_DNA_BSC1 = 'amber14/DNA.bsc1.xml'
    AMBER14_RNA_OL3 = 'amber14/RNA.OL3.xml'
    AMBER14_LIPID17 = 'amber14/lipid17.xml'
    AMBER14_GLYCAM_06J_1 = 'amber14/GLYCAM_06j-1.xml'
    CHARMM36 = 'charmm36.xml'


class WaterForceFields(Enum):
    AMBER14_TIP3P = 'amber14/tip3p.xml'
    AMBER14_TIP3PFB = 'amber14/tip3pfb.xml'
    AMBER14_TIP4PEW = 'amber14/tip4pew.xml'
    AMBER14_TIP4PFB = 'amber14/tip4pfb.xml'
    AMBER14_SPCE = 'amber14/spce.xml'
    AMBER14_OPC = 'amber14/opc.xml'
    AMBER14_OPC3 = 'amber14/opc3.xml'


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class GenGromacsFilesRequest(BaseModelMixin, BaseModel):
    pdbContent: str
    forceField: ForceFields
    waterForceField: WaterForceFields


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunSimulationsBase(BaseModelMixin, BaseModel):
    temperatureK: float = 273.15
    frictionCoeff: float = 1.0
    stepSize: float = 0.002
    integrator: Integrators = Integrators.langevin
    takeFrameEvery: int = 1000
    totalFrames: int = 10000


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunGromacsSimulationsRequest(RunSimulationsBase):
    top: str
    gro: str


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunPdbSimulationsRequest(RunSimulationsBase):
    pdbContent: str
    forceField: str
    waterForceField: WaterForceFields


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunPdbFixerRequest(BaseModelMixin, BaseModel):
    replaceNonStandardResidues: bool = False
    addMissingAtoms: bool = False
    addMissingHydrogens: bool = True
    pdbContent: str


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunPdbFixerResponse(BaseModelMixin, ErrorResponseMixing, BaseModel):
    pdbContent: str | None


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class RunSimulationsResponse(BaseModelMixin, ErrorResponseMixing, BaseModel):
    pdbContent: str | None


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class GenGroTopRequest(BaseModelMixin, BaseModel):
    ignoreMissingAtoms: bool = False
    forceField: ForceFields
    waterForceField: WaterForceFields
    pdbContent: str


@pydantic.dataclasses.dataclass
@dataclasses.dataclass(kw_only=True)
class GenGroTopResponse(BaseModelMixin, ErrorResponseMixing, BaseModel):
    gro: str | None
    top: str | None
