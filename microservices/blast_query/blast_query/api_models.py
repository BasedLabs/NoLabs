from pydantic import BaseModel

class SequenceQuery(BaseModel):
    """A query for a BLAST search.

    - sequence: could be nucleotide sequence for blastn, tblastx, or tblastn, or amino acid sequence for blastp or blastx.
    """
    sequence: str
    descriptions: int = 10
    alignments: int = 10
    hitlist_size: int = 10
    expect: float = 10.0
