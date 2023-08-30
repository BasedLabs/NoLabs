import os
from typing import List
from src.ai.model import BaseModel

class FakeFolding(BaseModel):

    def __init__(self, model_name: str, gpu: bool, model_task = ""):
        super().__init__(model_name, gpu, model_task)
        self.model_name = model_name
        self.gpu = gpu

    def load_model(self):
        pass

    def read_pdb_file(self, file_path):
        try:
            with open(file_path, 'r') as file:
                pdb_contents = file.read()
                return "".join(pdb_contents.split())
        except FileNotFoundError:
            print(f"Error: File '{file_path}' not found.")
            return None

    def predict(self, sequence: str) -> List[str]:
        dirname = os.path.dirname
        file_path = dirname(os.path.abspath(__file__)) + "/mock_resources/1r6a.pdb"
        mock_outputs = self.read_pdb_file(file_path)

        return mock_outputs