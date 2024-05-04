from __future__ import annotations
from typing import List, Optional
from uuid import UUID

from fastapi import UploadFile
from pydantic.dataclasses import dataclass


@dataclass
class RunProteinDesignRequest:
    job_id: UUID
    pdb_file: Optional[UploadFile] = None
    contig: str = '50'
    number_of_designs: int = 1
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None


@dataclass
class RunProteinDesignResponse:
    job_id: UUID
    pdb_files: List[str]


@dataclass
class JobPropertiesResponse:
    pdb_file: str
    pdb_file_name: str
    contig: str
    number_of_designs: int
    hotspots: Optional[str] = None
    timesteps: Optional[int] = None

    @staticmethod
    def default() -> JobPropertiesResponse:
        return JobPropertiesResponse(
            pdb_file='',
            pdb_file_name='',
            contig='',
            number_of_designs=2,
            hotspots=None,
            timesteps=None
        )


@dataclass
class GetJobResponse:
    job_id: UUID
    job_name: str
    pdb_files: List[str]
    properties: JobPropertiesResponse