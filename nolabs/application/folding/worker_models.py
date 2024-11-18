from pydantic import BaseModel


class InferenceInput(BaseModel):
    fasta_sequence: str

class InferenceOutput(BaseModel):
    pdb_content: str