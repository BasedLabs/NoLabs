from typing import Optional, List

from pydantic import BaseModel, Field


class RunRfdiffusionRequest(BaseModel):
    pdb_content: str
    contig: str
    hotspots: Optional[str] = ''
    inpaint: Optional[str] = ''
    timesteps: Optional[int] = 50
    number_of_designs: Optional[int] = 1
    remove_chain: Optional[str] = ''


class RunRfdiffusionResponse(BaseModel):
    pdbs_content: List[str] = Field(default_factory=list)
    errors: List[str] = Field(default_factory=list)