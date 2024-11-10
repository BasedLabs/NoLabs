from typing import Optional, List

from pydantic import BaseModel


class RunDiffDockPredictionInput(BaseModel):
    pdb_contents: str  # Protein contents in PDB format
    sdf_contents: str  # Ligand contents in SDF format
    inference_steps: int = 20
    samples_per_complex: int = 40
    batch_size: int = 10
    actual_steps: int = 18
    no_final_step_noise: bool = True

class SDFResult(BaseModel):
    sdf_file_name: str
    sdf_content: str
    scored_affinity: float
    minimized_affinity: float
    confidence: Optional[float] = None


class RunDiffDockPredictionOutput(BaseModel):
    pdb_contents: str
    sdf_results: List[SDFResult]