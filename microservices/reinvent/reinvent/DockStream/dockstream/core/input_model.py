from typing import List, Optional, Union

from pydantic import BaseModel

from dockstream.core.AutodockVina.AutodockVina_docker import AutodockVina
from dockstream.core.Corina.Corina_ligand_preparator import CorinaLigandPreparator
from dockstream.core.Gold.Gold_docker import Gold
from dockstream.core.OpenEye.OpenEye_docker import OpenEye
from dockstream.core.OpenEye.OpenEye_ligand_preparator import OpenEyeLigandPreparator
from dockstream.core.OpenEyeHybrid.OpenEyeHybrid_docker import OpenEyeHybrid
from dockstream.core.RDkit.RDkit_ligand_preparator import RDkitLigandPreparator
from dockstream.core.Schrodinger.Glide_docker import Glide
from dockstream.core.Schrodinger.Ligprep_ligand_preparator import LigprepLigandPreparator
from dockstream.core.rDock.rDock_docker import rDock


class EnvVariable(BaseModel):
    key: str
    value: str


class Environment(BaseModel):
    export: Optional[List[EnvVariable]]


class Logging(BaseModel):
    logfile: str


class Header(BaseModel):
    environment: Optional[Environment]
    logging: Logging


AnyLigandPreparator = Union[
    # CorinaLigandPreparator,
    LigprepLigandPreparator,
    # OpenEyeLigandPreparator,
    # RDkitLigandPreparator,
]


class LigandPreparation(BaseModel):
    """Ligand preparation: Specify embedding pool/ligand preparator."""
    embedding_pools: Union[AnyLigandPreparator, List[AnyLigandPreparator]]


AnyDocker = Union[
    # AutodockVina,
    Glide,
    # Gold,
    # OpenEye,
    # OpenEyeHybrid,
    # rDock,
]


class DockingInput(BaseModel):
    """Docking input.

    Consists of two big parts: ligand preparation and docking runs/backend.
    """

    header: Header
    ligand_preparation: LigandPreparation
    docking_runs: Union[AnyDocker, Optional[List[AnyDocker]]]


class AzdockInput(BaseModel):
    """Welcome to AZDock."""
    docking: DockingInput
