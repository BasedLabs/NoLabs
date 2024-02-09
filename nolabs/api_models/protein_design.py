import datetime
from typing import List, Optional

from fastapi import UploadFile
from pydantic import dataclasses


@dataclasses.dataclass
class RunProteinDesignRequest:
    experiment_name: str
    experiment_id: str
    pdb_file: Optional[UploadFile] = None
    contig: str = '50'
    number_of_designs: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None


@dataclasses.dataclass
class RunProteinDesignResponse:
    experiment_id: str
    experiment_name: str
    pdb_files: List[str]


@dataclasses.dataclass
class ExperimentPropertiesResponse:
    pdb_file: str
    pdb_file_name: str
    contig: str
    number_of_designs: int
    hotspots: Optional[str] = None
    timesteps: Optional[int] = None


@dataclasses.dataclass
class GetExperimentResponse:
    experiment_id: str
    experiment_name: str
    pdb_files: List[str]
    properties: ExperimentPropertiesResponse