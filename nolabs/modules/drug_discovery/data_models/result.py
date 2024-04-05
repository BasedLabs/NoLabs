from typing import List

from pydantic.dataclasses import dataclass

@dataclass
class JobId:
    value: str

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise ValueError()

        return self.value == other.value

    def __str__(self):
        return self.value

    def __hash__(self):
        return self.value.__hash__()

@dataclass
class UmolResultMetaData:
    job_id: str
    target_id: str
    ligand_id: str


@dataclass
class JobMetaData:
    job_id: str
    target_id: str
    ligand_id: str
    folding_method: str
    docking_method: str


@dataclass
class UmolDockingResultData:
    predicted_pdb: str
    predicted_sdf: str
    plddt_array: List[int]

@dataclass
class DiffDockDockingResultData:
    predicted_pdb: str
    predicted_sdf_file_name: str
    predicted_sdf_contents: str
    scored_affinity: float
    minimized_affinity: float
    confidence: float | None = None


@dataclass
class DiffDockResultMetaData:
    job_id: str
    target_id: str
    ligand_id: str
    ligand_file_name: str
    minimized_affinity: float
    scored_affinity: float
    confidence: float | None = None

@dataclass
class DiffDockLigandResultData:
    predicted_sdf_contents: str


@dataclass
class DiffDockModelParams:
    inference_steps: int = 20
    samples_per_complex: int = 1
    batch_size: int = 10
    actual_steps: int = 18