from pydantic import BaseModel


class PredictFoldingJobRequest(BaseModel):
    protein_sequence: str


class PredictFoldingJobResponse(BaseModel):
    pdb_content: str