from typing import Optional, List, Any

import pydantic
from fastapi import UploadFile
from pydantic import model_validator


@pydantic.dataclasses.dataclass
class RunAminoAcidRequest:
    experiment_name: str
    experiment_id: str
    amino_acid_sequence: Optional[str]
    fastas: Optional[List[UploadFile]]

    @classmethod
    @model_validator(mode='after')
    def check_inputs(cls, data: Any) -> Any:
        if not isinstance(data, RunAminoAcidRequest):
            raise ValueError('Incorrect data type')
        if not data.amino_acid_sequence and not data.fastas:
            raise ValueError('Either specify aminoacid sequence or fastas files')
        return data
