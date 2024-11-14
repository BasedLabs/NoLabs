from typing import Optional, List, Dict
from pydantic import BaseModel


class RunProteinMPNNPredictionInput(BaseModel):
    pdb_contents: str
    is_homomer: bool
    chains_to_design: Optional[List[str]] = None
    fixed_positions: Optional[Dict[str, List[int]]] = None
    num_seq_per_target: int = 2
    sampling_temp: float = 0.1
    seed: int = 37
    batch_size: int = 1


class SequenceResult(BaseModel):
    sequence: str
    fasta_content: str
    score: Optional[float] = None
    global_score: Optional[float] = None
    T: Optional[float] = None
    sample: Optional[int] = None
    seq_recovery: Optional[float] = None


class RunProteinMPNNPredictionOutput(BaseModel):
    sequences: List[SequenceResult]
