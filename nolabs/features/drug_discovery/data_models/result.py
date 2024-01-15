import dataclasses
from typing import List

from pydantic.dataclasses import dataclass

@dataclasses.dataclass
@dataclass
class ResultId:
    value: str

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise ValueError()

        return self.value == other.value

    def __str__(self):
        return self.value

    def __hash__(self):
        return self.value.__hash__()

@dataclasses.dataclass
@dataclass
class ResultMetaData:
    result_id: str
    target_id: str
    ligand_id: str

@dataclasses.dataclass
@dataclass
class DockingResultData:
    predicted_pdb: str
    predicted_sdf: str
    plddt_array: List[int]
