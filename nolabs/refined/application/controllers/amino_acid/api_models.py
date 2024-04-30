from typing import Optional, List, Any
from uuid import UUID

from fastapi import UploadFile
from pydantic import model_validator, BaseModel


class RunAminoAcidRequest(BaseModel):
    job_id: Optional[UUID]
    experiment_id: UUID
    fastas: List[UploadFile]

    @classmethod
    @model_validator(mode='after')
    def check_inputs(cls, data: Any) -> Any:
        if not data.fastas:
            raise ValueError('Specify at least one .fasta file')
        return data
