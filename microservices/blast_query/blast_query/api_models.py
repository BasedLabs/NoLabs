from enum import Enum
from typing import Optional

from pydantic import BaseModel


class BlastType(str, Enum):
    blastn = "blastn"
    blastp = "blastp"
    blastx = "blastx"
    tblastn = "tblastn"
    tblastx = "tblastx"


class SequenceQuery(BaseModel):
    """A query for a BLAST search.

    - sequence: could be nucleotide sequence for blastn, tblastx, or tblastn, or amino acid sequence for blastp or blastx.
    """

    sequence: str
    type: BlastType
    descriptions: Optional[int] = 10
    alignments: Optional[int] = 10
    hitlist_size: Optional[int] = 10
    expect: Optional[float] = 10.0
    job_id: Optional[str] = None


class IsJobRunningResponse(BaseModel):
    is_running: bool
