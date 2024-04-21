from typing import Optional, List, Any

from fastapi import UploadFile
from pydantic import model_validator


class RunAminoAcidRequest:
    experiment_name: str
    experiment_id: str
    fastas: Optional[List[UploadFile]]

    @classmethod
    @model_validator(mode='after')
    def check_inputs(cls, data: Any) -> Any:
        if not data.fastas:
            raise ValueError('Specify at least one .fasta file')
        return data
